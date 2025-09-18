/*  This is a program to remove a C-water watermark. */
/*  Copyright (C) 2025  Isaac Landry */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Affero General Public License as */
/* published by the Free Software Foundation, either version 3 of the */
/* License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but */
/* WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Affero General Public License for more details. */

/* You should have received a copy of the GNU Affero General Public */
/* License along with this program.  If not, see */
/* <https://www.gnu.org/licenses/>. */

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
verify_c_water(GimpDrawable* drawable,
                 guchar* pixels_to_change,
                 gint lower_limit_x,
                 gint lower_limit_y,
                 gint upper_limit_x,
                 gint upper_limit_y,
                 gint channels);

// global variables

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

  gimp_ui_init("verifyc_water", FALSE);
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

  gimp_progress_init("Verify C_Water...");

  gint x1, y1, x2, y2;

  gimp_drawable_mask_bounds(drawable->drawable_id, &x1, &y1, &x2, &y2);

  guchar* pixels_to_change;
  gint channels = gimp_drawable_bpp(drawable->drawable_id);
  pixels_to_change = g_new(guchar, channels * (x2 - x1) * (y2 - y1));

  verify_c_water(drawable, pixels_to_change, x1, x2, y1, y2, channels);

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

  gimp_install_procedure("plug-in-verify-c_water",
                         "verify-c_water",
                         "Displays \"Verify c_water\" in a dialog",
                         "Isaac Landry",
                         "Copyright Isaac Landry",
                         "2004",
                         "_Verify c_water...",
                         "RGB*, GRAY*",
                         GIMP_PLUGIN,
                         G_N_ELEMENTS(args),
                         0,
                         args,
                         NULL);

  gimp_plugin_menu_register("plug-in-verify-c_water", "<Image>/Filters/Misc");
}

static void
verify_c_water(GimpDrawable* drawable,
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
  guchar* c_water_bits; // The set of bits used for watermarking. We
                          // have 2 bits per 8 x 8 block.

  gimp_drawable_mask_bounds(drawable->drawable_id, &x1, &y1, &x2, &y2);
  width = x2 - x1;
  height = y2 - y1;

  gimp_pixel_rgn_init(&rgn_in, drawable, x1, y1, width, height, FALSE, FALSE);
  gimp_pixel_rgn_init(&rgn_out, drawable, x1, y1, width, height, TRUE, TRUE);

  const int total_cols = channels * width;

  for (int i = 0; i < 8; ++i) {
    row_arr[i] = g_new(guchar, total_cols);
  }

  int max_col_32 = total_cols - (total_cols % 32);
  printf("max_col_32: %d \n ", max_col_32);

  int max_row_8 = height - (height % 8);

  size_t num_blocks = (max_col_32 / 8) * (max_row_8 / 8);
  size_t c_water_bits_size = num_blocks / 4;
  // We have 3 channels (colors) per pixel, and 8 x 8 pixels per
  // block. We aim to get the number of blocks in the image. For 1024
  // x 1024 images, there are (1024 x 1024) * 3 / (8 x 8) = (2^14) * 3
  // blocks = 49152 blocks in this case.  Since there are 2 bits per
  // block, we have 49152 blocks * 2 bits / block = 98304 bits = 98304
  // bits / (8 bits / byte) = 12288 bytes.  The factor of division by
  // 4 results from (2 bits / block) / (8 bits / byte).  The size is
  // for byte addressing.
  c_water_bits = g_new(guchar, c_water_bits_size);

  guchar* edge_sign = g_new(guchar, num_blocks);
  guchar* central_sign = g_new(guchar, num_blocks);

  for (i = y1; i < y1 + max_row_8; i += 8) {
    /* Get row i through i+7 */
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[0], x1, i, width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[1], x1, MIN(y2 - 1, i + 1), width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[2], x1, MIN(y2 - 1, i + 2), width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[3], x1, MIN(y2 - 1, i + 3), width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[4], x1, MIN(y2 - 1, i + 4), width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[5], x1, MIN(y2 - 1, i + 5), width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[6], x1, MIN(y2 - 1, i + 6), width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[7], x1, MIN(y2 - 1, i + 7), width);

    // Break up into 8x8 subblocks of pixels

    for (gint col_offset = 0; col_offset < max_col_32; col_offset += 8) {
      int x_block = col_offset / 8; // the column index of each block.
      int y_block = (i - y1) / 8;   // the row index of each block
                                    // there are max_col_32 / 8 blocks per row.
      int block_index = x_block + y_block * max_col_32 / 8;

      // 2 bits per block are used. Gives the byte for the block.
      int c_water_bit_index = block_index / 4;

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

      int c_water_value =
        ((G_6_6_central >= 4 && G_6_6_central < 12) ? 1 : 0) +
        ((G_6_6_edge >= 4 && G_6_6_edge < 12) ? 2 : 0);

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
      c_water_bits[c_water_bit_index] =
        c_water_bits[c_water_bit_index] |
        (c_water_value << (sub_block_index * 2));
    }

    if (i % 10 == 0)
      gimp_progress_update((gdouble)(i - y1) / (gdouble)(height));
  }

  printf("blake3_hash verify: ");
  for (i = 0; i < BLAKE3_OUT_LEN; ++i) {
    printf("%.2x ", c_water_bits[i]);
  }
  printf("\n");

  // Check if the hash is too big to be in the c_water.
  if (BLAKE3_OUT_LEN > c_water_bits_size) {
    g_message("Could not recover the c_water because the hash is too big to be in the c_water");
    return;
  }

  unsigned char* bitmaps[1] = { c_water_bits };
  
  struct jbg_dec_state sd;

  jbg_dec_init(&sd);
  size_t dec_offset;

  int jbig_result = jbg_dec_in(&sd,
                               c_water_bits + BLAKE3_OUT_LEN,
                               c_water_bits_size - BLAKE3_OUT_LEN,
                               &dec_offset);

  printf("jbg_strerror: %s \n", jbg_strerror(jbig_result));
  printf("jbig_result: %d \n", jbig_result);

  int number_of_planes = jbg_dec_getplanes(&sd);

  // Recover the hash and the original bits from the compressed data.
  unsigned char* recovered_bits = jbg_dec_getimage(&sd, 0);

  unsigned long result_size = jbg_dec_getsize(&sd);
  printf("result_size: %ld \n", result_size);

  if (jbig_result != JBG_EOK) {
    g_message("No c_water detected in this image.\n");
    return;
  }

  unsigned long result_width = jbg_dec_getwidth(&sd);
  unsigned long result_height = jbg_dec_getheight(&sd);
  printf("result_width: %ld \n", result_width);
  printf("result_height: %ld \n", result_height);

  // Initialize the hasher.
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);

  // Restore the original image.
  for (i = y1; i < y1 + max_row_8; i += 8) {
    /* Get row i through i+7 */
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[0], x1, i, width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[1], x1, MIN(y2 - 1, i + 1), width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[2], x1, MIN(y2 - 1, i + 2), width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[3], x1, MIN(y2 - 1, i + 3), width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[4], x1, MIN(y2 - 1, i + 4), width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[5], x1, MIN(y2 - 1, i + 5), width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[6], x1, MIN(y2 - 1, i + 6), width);
    gimp_pixel_rgn_get_row(&rgn_in, row_arr[7], x1, MIN(y2 - 1, i + 7), width);

    // Break up into 8x8 subblocks of pixels
    for (gint col_offset = 0; col_offset < max_col_32; col_offset += 8) {
      int x_block = col_offset / 8; // the column index of each block.
      int y_block = (i - y1) / 8;   // the row index of each block
      int block_index =
        x_block +
        y_block * max_col_32 / 8; // there are max_col_32 / 8 blocks per row.

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

      int c_water_value =
        (c_water_bits[recovered_bit_index] >> (sub_block_index * 2)) & 3;
      // Select the new values of the bits for each pixel.
      //
      // new_value has two bits. "&" it with 1 to get the lowest
      // bit. divide by 2 to get the next bit.
      // The order is arbitrary.

      int c_water_bit_1_p_alpha = (c_water_value & 1);
      int c_water_bit_alpha = (c_water_value / 2);

      if (c_water_bit_1_p_alpha != bit_1_p_alpha) {
        // If we are adding a bit, then bit_1_p_alpha == 1 and
        // original_bit_1_p_alpha == 0.
        int add_subtract = (c_water_bit_1_p_alpha == 0) ? 1 : -1;
        add_subtract *= (central_sign[block_index] == 0) ? 1 : -1;
        for (x = 1; x <= 2; ++x) {
          for (y = 1; y <= 2; ++y) {
            int sgn = ((((x + y) % 2) == 0) ? 1 : -1) * add_subtract;

            row_arr[y][col_offset + x] =
              add_16(row_arr[y][col_offset + x], sgn);
          }
        }
      }

      if (c_water_bit_alpha != bit_alpha) {
        // If we are adding a bit, then bit_alpha == 1 and original_bit_alpha ==
        // 0.
        int add_subtract = (c_water_bit_alpha == 0) ? 1 : -1;
        add_subtract *= (edge_sign[block_index] == 0) ? 1 : -1;

        for (x = 0; x <= 3; x += 3) {
          for (y = 1; y <= 2; ++y) {
            int sgn = ((((x + y) % 2) == 0) ? 1 : -1) * add_subtract;

            row_arr[y][col_offset + x] =
              add_16(row_arr[y][col_offset + x], sgn);
          }
        }
      }
    }

    for (k = 0; k < 8; ++k) {
      gimp_pixel_rgn_set_row(&rgn_out, row_arr[k], x1, i + k, width);
      // Hash all of the pixels in the image.
      blake3_hasher_update(&hasher, row_arr[k], total_cols);
    }

    if (i % 10 == 0)
      gimp_progress_update((gdouble)(i - y1) / (gdouble)(height));
  }

  // Finalize the hash. BLAKE3_OUT_LEN is the default output length, 32 bytes.
  uint8_t blake3_hash[BLAKE3_OUT_LEN];
  blake3_hasher_finalize(&hasher, blake3_hash, BLAKE3_OUT_LEN);

  printf("blake3_hash: ");
  for (i = 0; i < BLAKE3_OUT_LEN; ++i) {
    printf("0x%.2x ", blake3_hash[i]);
  }
  printf("\n");

  for (i = 0; i < 8; ++i) {
    g_free(row_arr[i]);
  }

  gimp_drawable_flush(drawable);
  gimp_drawable_merge_shadow(drawable->drawable_id, TRUE);
  gimp_drawable_update(drawable->drawable_id, x1, y1, width, height);

  g_message("Image successfully restored.\n");
}
