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
// unsigned char compressed_data[10000]; // Compressed version of the bits used
// for watermarking.
// unsigned char* compressed_data;
size_t current_len = 0;

// Add pixel_value with the lowest bits mod 16.
char
add_16(char pixel_value, int adjustment)
{
  char high_bits = pixel_value & 0xF0;
  char low_bits = pixel_value & 0x0F;

  // Add and only take the low bits. This is equivalent to taking mod 16.
  low_bits += adjustment;
  low_bits = low_bits & 0x0F;

  // Gives the original high_bits and the new low_bits.
  return high_bits | low_bits;
}

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

  guchar* row_arr[8];     // The array of rows of pixels in an image.
  int offset = 128;       // To map the values from 0-255 to -128-127
  guchar* watermark_bits; // The set of bits used for watermarking. We
                          // have 2 bits per 8 x 8 block.

  gimp_drawable_mask_bounds(drawable->drawable_id, &x1, &y1, &x2, &y2);
  width = x2 - x1;
  height = y2 - y1;

  gimp_pixel_rgn_init(
    &rgn_in, drawable, x1, y1, x2 - x1, y2 - y1, FALSE, FALSE);
  gimp_pixel_rgn_init(&rgn_out, drawable, x1, y1, x2 - x1, y2 - y1, TRUE, TRUE);

  for (int i = 0; i < 8; ++i) {
    row_arr[i] = g_new(guchar, channels * (x2 - x1));
  }

  size_t num_blocks = channels * (width / 8) * (height / 8);
  size_t watermark_bits_size = num_blocks / 4;
  // We have 3 channels (colors) per pixel, and 8 x 8 pixels per
  // block. We aim to get the number of blocks in the image. For 1024
  // x 1024 images, there are (1024 x 1024) * 3 / (8 x 8) = (2^14) * 3
  // blocks = 49152 blocks in this case.  Since there are 2 bits per
  // block, we have 49152 blocks * 2 bits / block = 98304 bits = 98304
  // bits / (8 bits / byte) = 12288 bytes.  The factor of division by
  // 4 results from (2 bits / block) / (8 bits / byte).  The size is
  // for byte addressing.
  watermark_bits = g_new(guchar, watermark_bits_size);

  guchar* edge_sign = g_new(guchar, num_blocks);
  guchar* central_sign = g_new(guchar, num_blocks);

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
      int x_block = col_offset / 8; // the column index of each block.
      int y_block =
        (i - y1) / 8; // the row index of each block
                      // there are channel * width / 8 blocks per row.
      int block_index = x_block + y_block * channels * width / 8;

      // 2 bits per block are used. Gives the byte for the block.
      int watermark_bit_index = block_index / 4;

      // Gives the bits within the byte for this block.
      int sub_block_index = block_index & 3;

      int32_t G_6_6_central = 0;

      for (x = 1; x <= 2; ++x) {
        for (y = 1; y <= 2; ++y) {
          int sgn = (((x + y) % 2) == 0) ? 1 : -1;
          G_6_6_central += sgn * (row_arr[y][col_offset + x] - offset);
          G_6_6_central += -sgn * (row_arr[y + 4][col_offset + x] - offset);
          G_6_6_central += -sgn * (row_arr[y][col_offset + x + 4] - offset);
          G_6_6_central += sgn * (row_arr[y + 4][col_offset + x + 4] - offset);
        }
      }

      int32_t G_6_6_edge = 0;

      for (x = 0; x <= 3; x += 3) {
        for (y = 1; y <= 2; ++y) {
          int sgn = (((x + y) % 2) == 0) ? 1 : -1;
          G_6_6_edge += sgn * (row_arr[y][col_offset + x] +
                               row_arr[x][col_offset + y] - 2 * offset);
          G_6_6_edge += -sgn * (row_arr[y + 4][col_offset + x] +
                                row_arr[x][col_offset + y + 4] - 2 * offset);
          G_6_6_edge += -sgn * (row_arr[y][col_offset + x + 4] +
                                row_arr[x + 4][col_offset + y] - 2 * offset);
          G_6_6_edge += sgn * (row_arr[y + 4][col_offset + x + 4] +
                               row_arr[x + 4][col_offset + y + 4] - 2 * offset);
        }
      }

      // Pack the bits to save into orig_value.
      G_6_6_central = (G_6_6_central + 8192) % 16;
      G_6_6_edge = (G_6_6_edge + 8192) % 16;

      central_sign[block_index] = (G_6_6_central >= 8) ? 1 : 0;
      edge_sign[block_index] = (G_6_6_edge >= 8) ? 1 : 0;

      int watermark_value =
        ((G_6_6_central >= 4 && G_6_6_central < 12) ? 1 : 0) +
        ((G_6_6_edge >= 4 && G_6_6_edge < 12) ? 2 : 0);

      if (x_block < 16 && y_block == 0) {
        printf("G: %d %d %d %d %d %d \n",
               x_block,
               y_block,
               block_index,
               G_6_6_central,
               G_6_6_edge,
               watermark_value);
      }

      // Save the bits into the original_bits array.
      //
      // Get the byte to set with original_bits[original_bit_index].
      //
      // Shift orig_value using "(orig_value << (sub_block_index * 2))"
      // so that the two bits to set are in the right place in the byte.
      // (we multiply sub_block_index by 2 because there are 2 bits
      // per block).
      //
      // We "|" these two together to update the bit pattern for
      // that byte.
      watermark_bits[watermark_bit_index] =
        watermark_bits[watermark_bit_index] |
        (watermark_value << (sub_block_index * 2));
    }

    if (i % 10 == 0)
      gimp_progress_update((gdouble)(i - y1) / (gdouble)(y2 - y1));
  }

  printf("blake3_hash verify: ");
  for (i = 0; i < BLAKE3_OUT_LEN; ++i) {
    printf("%.2x ", watermark_bits[i]);
  }
  printf("\n");

  printf("jbg: ");
  for (i = BLAKE3_OUT_LEN; i < BLAKE3_OUT_LEN + 32; ++i) {
    printf("%.2x ", watermark_bits[i]);
  }
  printf("\n");

  unsigned char* bitmaps[1] = { watermark_bits };
  // TODO: this shouldn't be bigger than 2 * watermark_bits_size. Otherwise, it
  // won't fit. compressed_data = malloc(2 * watermark_bits_size);

  struct jbg_dec_state sd;

  jbg_dec_init(&sd);
  size_t dec_offset;
  int jbig_result = jbg_dec_in(&sd,
                               watermark_bits + BLAKE3_OUT_LEN,
                               watermark_bits_size - BLAKE3_OUT_LEN,
                               &dec_offset);

  printf("offset: %d \n", dec_offset);
  printf("jbig_result: %d \n", jbig_result);
  printf("JBG_EAGAIN: %d \n", JBG_EAGAIN);
  printf("JBG_EOK: %d \n", JBG_EOK);
  printf("JBG_EOK_INTR: %d \n", JBG_EOK_INTR);

  int number_of_planes = jbg_dec_getplanes(&sd);

  // Recover the hash and the original bits from the compressed data.
  unsigned char* recovered_bits = jbg_dec_getimage(&sd, 0);
  long result_size = jbg_dec_getsize(&sd);

  printf("Result size: %d \n", result_size);
  printf("Width: %d \n", jbg_dec_getwidth(&sd));
  printf("Height: %d \n", jbg_dec_getheight(&sd));

  printf("original_bit_size: %d \n", watermark_bits_size);

  // Initialize the hasher.
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);

  // Restore the original image.
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
      int x_block = col_offset / 8; // the column index of each block.
      int y_block = (i - y1) / 8;   // the row index of each block
      int block_index =
        x_block + y_block * channels * width /
                    8; // there are channel * width / 8 blocks per row.

      // 2 bits per block are used. Gives the byte for the block.
      int recovered_bit_index = block_index / 4;

      // Gives the bits within the byte for this block.
      int sub_block_index = block_index & 3;

      // 1_p_alpha means 1 + alpha
      // TODO: Hardcode which pixels to change. This should be fixed
      // to use the hash to choose which pixels.
      const int x_bit_1_p_alpha_index = 0;
      const int y_bit_1_p_alpha_index = 0;

      const int x_bit_alpha_index = 1;
      const int y_bit_alpha_index = 0;

      // Get the bits that we'll be setting.
      //
      // Get the byte with new_bits[original_bit_index].
      //
      // Shift the byte so that the two bits we want are all the way
      // to the right with ">> (sub_block_index * 2))" (we multiply
      // sub_block_index by 2 because there are 2 bits per block).
      //
      // We use "& 3" to select the lowest two bits because other
      // bits are for other blocks.
      int recovered_value =
        (recovered_bits[recovered_bit_index] >> (sub_block_index * 2)) & 3;

      // Select the new values of the bits for each pixel.
      //
      // new_value has two bits. "&" it with 1 to get the lowest
      // bit. divide by 2 to get the next bit.
      // The order is arbitrary.
      int bit_1_p_alpha = (recovered_value & 1);
      int bit_alpha = (recovered_value / 2);

      int watermark_value =
        (watermark_bits[recovered_bit_index] >> (sub_block_index * 2)) & 3;
      // Select the new values of the bits for each pixel.
      //
      // new_value has two bits. "&" it with 1 to get the lowest
      // bit. divide by 2 to get the next bit.
      // The order is arbitrary.

      int watermark_bit_1_p_alpha = (watermark_value & 1);
      int watermark_bit_alpha = (watermark_value / 2);

      /* // The 1_p_alpha pixel with the lowest bit set to 0. */
      /* int row_1_p_alpha = */
      /*   (row_arr[y_bit_1_p_alpha_index][x_bit_1_p_alpha_index + col_offset] &
       */
      /*    0xfe); */

      /* // The alpha pixel with the lowest bit set to 0. */
      /* int row_alpha = */
      /*   (row_arr[y_bit_alpha_index][x_bit_alpha_index + col_offset] & 0xfe);
       */

      if (watermark_bit_1_p_alpha != bit_1_p_alpha) {
        // If we are adding a bit, then bit_1_p_alpha == 1 and
        // original_bit_1_p_alpha == 0.
        int add_subtract = (watermark_bit_1_p_alpha == 0) ? 1 : -1;
        add_subtract *= (central_sign[block_index] == 0) ? 1 : -1;
        for (x = 1; x <= 2; ++x) {
          for (y = 1; y <= 2; ++y) {
            int sgn = ((((x + y) % 2) == 0) ? 1 : -1) * add_subtract;

            /* if (x_block >= 72 && x_block < 90 && y_block == 0) { */
            /*   printf("1_p_alpha: %d %d %d %d \n", */
            /*          x, */
            /*          y, */
            /*          sgn, */
            /*          row_arr[y][col_offset + x]); */
            /* } */

            row_arr[y][col_offset + x] =
              add_16(row_arr[y][col_offset + x], sgn);

            /* row_arr[y + 4][col_offset + x] = */
            /*   add_16(row_arr[y + 4][col_offset + x], -sgn); */
            /* row_arr[y][col_offset + x + 4] = */
            /*   add_16(row_arr[y][col_offset + x + 4], -sgn); */
            /* row_arr[y + 4][col_offset + x + 4] = */
            /*   add_16(row_arr[y + 4][col_offset + x + 4], sgn); */
          }
        }
      }

      if (watermark_bit_alpha != bit_alpha) {
        // If we are adding a bit, then bit_alpha == 1 and original_bit_alpha ==
        // 0.
        int add_subtract = (watermark_bit_alpha == 0) ? 1 : -1;
        add_subtract *= (edge_sign[block_index] == 0) ? 1 : -1;

        for (x = 0; x <= 3; x += 3) {
          for (y = 1; y <= 2; ++y) {
            int sgn = ((((x + y) % 2) == 0) ? 1 : -1) * add_subtract;

            if (x_block == 0 && y_block == 0) {
              printf("alpha: %d %d %d %d %d %d %d\n",
                     x,
                     y,
                     sgn,
                     add_subtract,
                     watermark_bit_alpha,
                     edge_sign[block_index],
                     row_arr[y][col_offset + x]);
            }

            row_arr[y][col_offset + x] =
              add_16(row_arr[y][col_offset + x], sgn);
            /* row_arr[y + 4][col_offset + x] = */
            /*   add_16(row_arr[y + 4][col_offset + x], -sgn); */
            /* row_arr[y][col_offset + x + 4] = */
            /*   add_16(row_arr[y][col_offset + x + 4], -sgn); */
            /* row_arr[y + 4][col_offset + x + 4] = */
            /*   add_16(row_arr[y + 4][col_offset + x + 4], sgn); */

            /* row_arr[x][col_offset + y] = */
            /*   add_16(row_arr[x][col_offset + y], sgn); */
            /* row_arr[x + 4][col_offset + y] = */
            /*   add_16(row_arr[x + 4][col_offset + y], -sgn); */
            /* row_arr[x][col_offset + y + 4] = */
            /*   add_16(row_arr[x][col_offset + y + 4], -sgn); */
            /* row_arr[x + 4][col_offset + y + 4] = */
            /*   add_16(row_arr[x + 4][col_offset + y + 4], sgn); */
          }
        }
      }

      if (i == y1 + 128 && col_offset == 8) {
        printf("\nbit_1_p_alpha: %d \n", bit_1_p_alpha);
        printf("bit_alpha: %d \n", bit_alpha);
        printf("watermark_bit_1_p_alpha: %d \n", watermark_bit_1_p_alpha);
        printf("watermark_bit_alpha: %d \n", watermark_bit_alpha);
        printf("recovered_value: %d \n", recovered_value);
        printf("recovered_bit_index: %d \n", recovered_bit_index);
        printf("recovered_bits: %.2x \n", recovered_bits[recovered_bit_index]);

        printf(
          "Before: %.2x %.2x\n",
          row_arr[y_bit_1_p_alpha_index][x_bit_1_p_alpha_index + col_offset],
          row_arr[y_bit_alpha_index][x_bit_alpha_index + col_offset]);
      }

      /* // Add the new bit value to the 1_p_alpha pixel. */
      /* row_arr[y_bit_1_p_alpha_index][x_bit_1_p_alpha_index + col_offset] = */
      /*   (row_1_p_alpha + bit_1_p_alpha); */

      /* // Add the new bit value to the alpha pixel. */
      /* row_arr[y_bit_alpha_index][x_bit_alpha_index + col_offset] = */
      /*   (row_alpha + bit_alpha); */

      if (i == y1 + 128 && col_offset == 8) {
        printf(
          "After: %.2x %.2x\n",
          row_arr[y_bit_1_p_alpha_index][x_bit_1_p_alpha_index + col_offset],
          row_arr[y_bit_alpha_index][x_bit_alpha_index + col_offset]);
      }
    }

    for (k = 0; k < 8; ++k) {
      gimp_pixel_rgn_set_row(&rgn_out, row_arr[k], x1, i + k, x2 - x1);
      // Hash all of the pixels in the image.
      blake3_hasher_update(&hasher, row_arr[k], channels * (x2 - x1));
    }

    /* if (i <= y1 + 127) { */
    /*   /\* printf("first row of row_arr: "); *\/ */
    /*   /\* for (j = x1; (j < x2) && (j < x1 + 100); ++j){ *\/ */
    /*   /\* 	printf("%.2x ", (row_arr[0][j])); *\/ */
    /*   /\* } *\/ */
    /*   /\* printf("\n"); *\/ */

    /*   for (k = 0; k < 8; ++k) { */
    /*     // Hash all of the pixels in the image. */
    /*     blake3_hasher_update(&hasher, row_arr[k], channels * (x2 - x1)); */
    /*   } */
    /* } */

    if (i == y1 + 128) {
      printf("128'th row of row_arr: ");
      for (j = x1; (j < x2) && (j < x1 + 100); ++j) {
        printf("%.2x ", (row_arr[0][j]));
      }
      printf("\n");
    }

    if (i % 10 == 0)
      gimp_progress_update((gdouble)(i - y1) / (gdouble)(y2 - y1));
  }

  // Finalize the hash. BLAKE3_OUT_LEN is the default output length, 32 bytes.
  uint8_t blake3_hash[BLAKE3_OUT_LEN];
  blake3_hasher_finalize(&hasher, blake3_hash, BLAKE3_OUT_LEN);

  printf("blake3_hash: ");
  for (i = 0; i < BLAKE3_OUT_LEN; ++i) {
    printf("%.2x ", blake3_hash[i]);
  }
  printf("\n");

  for (i = 0; i < 8; ++i) {
    g_free(row_arr[i]);
  }

  printf("Assertion point.\n");

  gimp_drawable_flush(drawable);
  gimp_drawable_merge_shadow(drawable->drawable_id, TRUE);
  gimp_drawable_update(drawable->drawable_id, x1, y1, x2 - x1, y2 - y1);
}
