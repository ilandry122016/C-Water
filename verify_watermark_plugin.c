#include "blake3.h"
#include <errno.h>
#include <jbig.h>
#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

static void
query(void);
static void
run(const gchar* name,
    gint nparams,
    const GimpParam* param,
    gint* nreturn_vals,
    GimpParam** return_vals);

void
gimp_pixel_rgn_get_row(GimpPixelRgn* pr,
                       guchar* buf,
                       gint x,
                       gint y,
                       gint width);
void
gimp_pixel_rgn_set_row(GimpPixelRgn* pr,
                       const guchar* buf,
                       gint x,
                       gint y,
                       gint width);

static void
verify_watermark(GimpDrawable* drawable,
                 guchar* pixels_to_change,
                 gint lower_limit_x,
                 gint lower_limit_y,
                 gint upper_limit_x,
                 gint upper_limit_y,
                 gint channels);

// global variables
unsigned char compressed_data[10000];
size_t current_len = 0;

GimpPlugInInfo PLUG_IN_INFO = { NULL, NULL, query, run };

MAIN()

static void
run(const gchar* name,
    gint nparams,
    const GimpParam* param,
    gint* nreturn_vals,
    GimpParam** return_vals)
{
  static GimpParam values[1];
  GimpPDBStatusType status = GIMP_PDB_SUCCESS;
  GimpRunMode run_mode;
  GimpDrawable* drawable;

  gimp_ui_init("verifywatermark", FALSE);
  GtkWidget* main_vbox;
  main_vbox = gtk_vbox_new(FALSE, 6);
  gtk_widget_show(main_vbox);

  /* Setting mandatory output values */
  *nreturn_vals = 1;
  *return_vals = values;

  values[0].type = GIMP_PDB_STATUS;
  values[0].data.d_status = status;

  /* Getting run_mode - we won't display a dialog if
   * we are in NONINTERACTIVE mode */
  run_mode = param[0].data.d_int32;

  /*  Get the specified drawable  */
  drawable = gimp_drawable_get(param[2].data.d_drawable);

  gimp_progress_init("Verify Watermark...");

  gint x1, y1, x2, y2;

  gimp_drawable_mask_bounds(drawable->drawable_id, &x1, &y1, &x2, &y2);

  guchar* pixels_to_change;
  gint channels = gimp_drawable_bpp(drawable->drawable_id);
  pixels_to_change = g_new(guchar, channels * (x2 - x1) * (y2 - y1));

  verify_watermark(drawable, pixels_to_change, x1, x2, y1, y2, channels);

  gimp_displays_flush();
  gimp_drawable_detach(drawable);

  return;
}

static void
query(void)
{
  static GimpParamDef args[] = {
    { GIMP_PDB_INT32, "run-mode", "Run mode" },
    { GIMP_PDB_IMAGE, "image", "Input image" },
    { GIMP_PDB_DRAWABLE, "drawable", "Input drawable" }
  };

  gimp_install_procedure("plug-in-verify-watermark",
                         "verify-watermark",
                         "Displays \"Verify watermark\" in a dialog",
                         "Isaac Landry",
                         "Copyright Isaac Landry",
                         "2004",
                         "_Verify watermark...",
                         "RGB*, GRAY*",
                         GIMP_PLUGIN,
                         G_N_ELEMENTS(args),
                         0,
                         args,
                         NULL);

  gimp_plugin_menu_register("plug-in-verify-watermark", "<Image>/Filters/Misc");
}

static void
verify_watermark(GimpDrawable* drawable,
                 guchar* pixels_to_change,
                 gint lower_limit_x,
                 gint lower_limit_y,
                 gint upper_limit_x,
                 gint upper_limit_y,
                 gint channels)
{
  gint i, j, k;

  gint x1 = lower_limit_x;
  gint y1 = lower_limit_y;
  gint x2 = upper_limit_x;
  gint y2 = upper_limit_y;

  GimpPixelRgn rgn_in, rgn_out;
  gint width, height;

  gint u, v, x, y;

  int encode_u = 6;
  int encode_v = 6;

  guchar* row_arr[8];
  int offset = 128; // To map the values from 0-255 to -128-127
  guchar* original_bits;

  gimp_drawable_mask_bounds(drawable->drawable_id, &x1, &y1, &x2, &y2);
  width = x2 - x1;
  height = y2 - y1;

  gimp_pixel_rgn_init(
    &rgn_in, drawable, x1, y1, x2 - x1, y2 - y1, FALSE, FALSE);
  gimp_pixel_rgn_init(&rgn_out, drawable, x1, y1, x2 - x1, y2 - y1, TRUE, TRUE);

  for (int i = 0; i < 8; ++i) {
    row_arr[i] = g_new(guchar, channels * (x2 - x1));
  }

  size_t original_bits_size = channels * (width / 8) * (height / 8) / 4;
  original_bits = g_new(
    guchar, original_bits_size); // We need 2 bits per 8x8 block. One is for the
                                 // signal, and the other is for the modulator.

  for (i = y1; i < y2; i += 8) {
    /* Get row i through i+7 */
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[0], x1, i, x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[1], x1, MIN(y2 - 1, i + 1), x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[2], x1, MIN(y2 - 1, i + 2), x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[3], x1, MIN(y2 - 1, i + 3), x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[4], x1, MIN(y2 - 1, i + 4), x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[5], x1, MIN(y2 - 1, i + 5), x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[6], x1, MIN(y2 - 1, i + 6), x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[7], x1, MIN(y2 - 1, i + 7), x2 - x1);

    // Break up into 8x8 subblocks of pixels

    for (gint col_offset = 0; col_offset < channels * (x2 - x1);
         col_offset += 8) {
      int x_block = col_offset / 8;
      int y_block = (i - y1) / 8;
      int block_index = x_block + y_block * channels * width / 8;
      int original_bit_index = block_index / 4;
      int sub_block_index = block_index & 3;

      const int x_bit_1_p_alpha_index = 0;
      const int y_bit_1_p_alpha_index = 0;

      const int x_bit_alpha_index = 1;
      const int y_bit_alpha_index = 0;

      // Get the lowest 2 integer bits of the pixels.
      int bit_1_p_alpha =
        (row_arr[y_bit_1_p_alpha_index][x_bit_1_p_alpha_index + col_offset] &
         1);
      int bit_alpha =
        (row_arr[y_bit_alpha_index][x_bit_alpha_index + col_offset] & 1);

      int orig_value = bit_1_p_alpha + 2 * bit_alpha;
      original_bits[original_bit_index] = original_bits[original_bit_index] |
                                          (orig_value << (sub_block_index * 2));
    }

    if (i % 10 == 0)
      gimp_progress_update((gdouble)(i - y1) / (gdouble)(y2 - y1));
  }

  unsigned char* bitmaps[1] = { original_bits };

  struct jbg_dec_state sd;

  jbg_dec_init(&sd);
  size_t dec_offset;
  jbg_dec_in(&sd,
             original_bits + BLAKE3_OUT_LEN,
             original_bits_size - BLAKE3_OUT_LEN,
             &dec_offset);
  printf("offset: %d \n", dec_offset);

  printf("crash point 1\n");
  int number_of_planes = jbg_dec_getplanes(&sd);
  printf("crash point 2\n");
  unsigned char* result_bitmap = jbg_dec_getimage(&sd, 0);
  printf("crash point 3\n");
  long result_size = jbg_dec_getsize(&sd);

  printf("Result size: %d \n", result_size);
  printf("Width: %d \n", jbg_dec_getwidth(&sd));
  printf("Height: %d \n", jbg_dec_getheight(&sd));

  printf("original_bit_size: %d \n", original_bits_size);

  printf("blake3_hash: ");
  for (i = 0; i < BLAKE3_OUT_LEN; ++i) {
    printf("%.2x ", original_bits[i]);
  }
  printf("\n");

  for (i = y1; i < y2; i += 8) {
    /* Get row i through i+7 */
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[0], x1, i, x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[1], x1, MIN(y2 - 1, i + 1), x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[2], x1, MIN(y2 - 1, i + 2), x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[3], x1, MIN(y2 - 1, i + 3), x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[4], x1, MIN(y2 - 1, i + 4), x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[5], x1, MIN(y2 - 1, i + 5), x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[6], x1, MIN(y2 - 1, i + 6), x2 - x1);
    gimp_pixel_rgn_get_row(
      &rgn_in, row_arr[7], x1, MIN(y2 - 1, i + 7), x2 - x1);

    // Break up into 8x8 subblocks of pixels
    for (gint col_offset = 0; col_offset < channels * (x2 - x1);
         col_offset += 8) {
      int x_block = col_offset / 8;
      int y_block = (i - y1) / 8;
      int block_index = x_block + y_block * channels * width / 8;
      int original_bit_index = block_index / 4;
      int sub_block_index = block_index & 3;

      // int new_value = (new_bits[original_bit_index] >> (sub_block_index * 2))
      // & 3;
      int new_value = 0;
    }

    for (k = 0; k < 8; ++k) {
      gimp_pixel_rgn_set_row(&rgn_out, row_arr[k], x1, i + k, x2 - x1);
    }

    if (i % 10 == 0)
      gimp_progress_update((gdouble)(i - y1) / (gdouble)(y2 - y1));
  }

  for (i = 0; i < 8; ++i) {
    g_free(row_arr[i]);
  }

  /* gimp_pixel_rgn_set_row(&rgn_out, */
  /* 			 original_bits[1], */
  /* 			 x1, 1, */
  /* 			 x2 - x1); */

  printf("Assertion point.\n");

  gimp_drawable_flush(drawable);
  gimp_drawable_merge_shadow(drawable->drawable_id, TRUE);
  gimp_drawable_update(drawable->drawable_id, x1, y1, x2 - x1, y2 - y1);
}
