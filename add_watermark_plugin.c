#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include <jbig.h>
#include "blake3.h"
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>

static void query (void);
static void run (const gchar *name,
		 gint nparams,
		 const GimpParam *param,
		 gint *nreturn_vals,
		 GimpParam **return_vals);

void gimp_pixel_rgn_get_row (GimpPixelRgn *pr,
			     guchar *buf,
			     gint x,
			     gint y,
			     gint width);
void gimp_pixel_rgn_set_row (GimpPixelRgn *pr,
			     const guchar *buf,
			     gint x,
			     gint y,
			     gint width);

static void add_watermark (GimpDrawable *drawable, guchar *pixels_to_change, gint lower_limit_x, gint lower_limit_y, gint upper_limit_x, gint upper_limit_y, gint channels);

// global variables
unsigned char compressed_data[10000]; // Compressed version of the bits used for watermarking.
size_t current_len = 0;

void output_bie(unsigned char *start, size_t len, void *file)
{
  for (int i = 0; i < len; ++i){
    compressed_data[current_len + i] = start[i];
  }
  current_len += len;
  
  return;
}

GimpPlugInInfo PLUG_IN_INFO = {
  NULL,
  NULL,
  query,
  run
};

MAIN()

  static void
  run (const gchar      *name,
       gint              nparams,
       const GimpParam  *param,
       gint             *nreturn_vals,
       GimpParam       **return_vals)
{
  static GimpParam  values[1];
  GimpPDBStatusType status = GIMP_PDB_SUCCESS;
  GimpRunMode       run_mode;
  GimpDrawable     *drawable;
	
  gimp_ui_init("addwatermark", FALSE);
  GtkWidget *main_vbox;
  main_vbox = gtk_vbox_new (FALSE, 6);
  gtk_widget_show (main_vbox);

  /* Setting mandatory output values */
  *nreturn_vals = 1;
  *return_vals  = values;

  values[0].type = GIMP_PDB_STATUS;
  values[0].data.d_status = status;

  /* Getting run_mode - we won't display a dialog if 
   * we are in NONINTERACTIVE mode */
  run_mode = param[0].data.d_int32;

  /*  Get the specified drawable  */
  drawable = gimp_drawable_get (param[2].data.d_drawable);
	
  gimp_progress_init ("Add Watermark...");

  gint x1, y1, x2, y2;
	
  gimp_drawable_mask_bounds (drawable->drawable_id,
                             &x1, &y1,
                             &x2, &y2);

  guchar *pixels_to_change;
  gint channels = gimp_drawable_bpp (drawable->drawable_id);
  pixels_to_change = g_new(guchar, channels * (x2 - x1) * (y2 - y1));

  add_watermark (drawable, pixels_to_change, x1, x2, y1, y2, channels);

  gimp_displays_flush ();
  gimp_drawable_detach (drawable);
	
  return;
}

static void
query (void)
{
  static GimpParamDef args[] = {
    {
      GIMP_PDB_INT32,
      "run-mode",
      "Run mode"
    },
    {
      GIMP_PDB_IMAGE,
      "image",
      "Input image"
    },
    {
      GIMP_PDB_DRAWABLE,
      "drawable",
      "Input drawable"
    }
  };

  gimp_install_procedure (
			  "plug-in-add-watermark",
			  "add-watermark",
			  "Displays \"Add watermark\" in a dialog",
			  "Isaac Landry",
			  "Copyright Isaac Landry",
			  "2004",
			  "_Add watermark...",
			  "RGB*, GRAY*",
			  GIMP_PLUGIN,
			  G_N_ELEMENTS (args), 0,
			  args, NULL);

  gimp_plugin_menu_register ("plug-in-add-watermark",
			     "<Image>/Filters/Misc");
}

static void
add_watermark(GimpDrawable *drawable, guchar *pixels_to_change, gint lower_limit_x, gint lower_limit_y, gint upper_limit_x, gint upper_limit_y, gint channels)
{
  gint         i, j, k;
  
  gint x1 = lower_limit_x;
  gint y1 = lower_limit_y;
  gint x2 = upper_limit_x;
  gint y2 = upper_limit_y;
  
  GimpPixelRgn rgn_in, rgn_out;
  gint         width, height;

  gint u, v, x, y;

  int encode_u = 6;
  int encode_v = 6;

  guchar* row_arr[8]; // The array of rows of pixels in an image.
  int offset = 128; // To map the values from 0-255 to -128-127
  guchar* original_bits; // The set of bits used for watermarking. We
			 // have 2 bits per 8 x 8 block.
  
  gimp_drawable_mask_bounds (drawable->drawable_id,
                             &x1, &y1,
                             &x2, &y2);
  width = x2 - x1;
  height = y2 - y1;
    
  gimp_pixel_rgn_init (&rgn_in,
                       drawable,
                       x1, y1,
                       x2 - x1, y2 - y1, 
                       FALSE, FALSE);
  gimp_pixel_rgn_init (&rgn_out,
                       drawable,
                       x1, y1,
                       x2 - x1, y2 - y1,
                       TRUE, TRUE);

  for (int i = 0; i < 8; ++i){
    row_arr[i] = g_new(guchar, channels * (x2 - x1));
  }

  size_t original_bits_size = channels * (width / 8) * (height / 8 ) / 4;
  // We have 3 channels (colors) per pixel, and 8 x 8 pixels per
  // block. We aim to get the number of blocks in the image. For 1024
  // x 1024 images, there are (1024 x 1024) * 3 / (8 x 8) = (2^14) * 3
  // blocks = 49152 blocks in this case.  Since there are 2 bits per
  // block, we have 49152 blocks * 2 bits / block = 98304 bits = 98304
  // bits / (8 bits / byte) = 12288 bytes.  The factor of division by
  // 4 results from (2 bits / block) / (8 bits / byte).  The size is
  // for byte addressing.
  original_bits = g_new(guchar, original_bits_size);

  // Initialize the hasher.
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);

  for (i = y1; i < y2; i += 8)
    {
      /* Get row i through i+7 */
      gimp_pixel_rgn_get_row (&rgn_in,
                              row_arr[0],
                              x1, i,
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
                              row_arr[1],
                              x1, MIN (y2 - 1, i + 1),
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
			      row_arr[2],
                              x1, MIN (y2 - 1, i + 2),
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
			      row_arr[3],
                              x1, MIN (y2 - 1, i + 3),
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
			      row_arr[4],
                              x1, MIN (y2 - 1, i + 4),
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
			      row_arr[5],
                              x1, MIN (y2 - 1, i + 5),
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
			      row_arr[6],
                              x1, MIN (y2 - 1, i + 6),
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
                              row_arr[7],
                              x1, MIN (y2 - 1, i + 7),
                              x2 - x1);

      for (j = 0; j < 8; ++j){
	blake3_hasher_update(&hasher, row_arr[j], channels * (x2 - x1));
	// Hash all of the pixels in the image.
      }

      // Break up into 8x8 subblocks of pixels
      
      for (gint col_offset = 0; col_offset < channels * (x2 - x1); col_offset += 8){
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
	int bit_1_p_alpha = (row_arr[y_bit_1_p_alpha_index][x_bit_1_p_alpha_index + col_offset] & 1);
	int bit_alpha = (row_arr[y_bit_alpha_index][x_bit_alpha_index + col_offset] & 1);
      
	int orig_value = bit_1_p_alpha + 2 * bit_alpha;
	original_bits[original_bit_index] = original_bits[original_bit_index] | (orig_value << (sub_block_index * 2));
      }      

      if (i % 10 == 0)
	gimp_progress_update ((gdouble) (i - y1) / (gdouble) (y2 - y1));
    }

  // Finalize the hash. BLAKE3_OUT_LEN is the default output length, 32 bytes.
  uint8_t blake3_hash[BLAKE3_OUT_LEN];
  blake3_hasher_finalize(&hasher, blake3_hash, BLAKE3_OUT_LEN);

  printf("blake3_hash: ");
  for (i = 0; i < BLAKE3_OUT_LEN; ++i){
    printf("%.2x ", blake3_hash[i]);
  }
  printf("\n");

  unsigned char *bitmaps[1] = {original_bits};
  struct jbg_enc_state se;
 
  jbg_enc_init(&se, (channels * width / 8) / 4, height / 8, 1, bitmaps, 
	       output_bie, stdout);              /* initialize encoder */
  jbg_enc_out(&se);                                    /* encode image */
  jbg_enc_free(&se);                    /* release allocated resources */

  printf("current_len: %d \n", current_len);
  printf("original_bits_size: %d \n", original_bits_size);
  
  guchar* new_bits = g_new(guchar, original_bits_size);
  // guchar is maybe unsigned char?

  // The watermark is the hash of all pixels in the image and the bits
  // before modification used for watermarking in a compressed form.

  // Store the hash of all pixels in the image.
  memcpy(new_bits, blake3_hash, BLAKE3_OUT_LEN);

  // Store the bits used for watermarking (2 pixels for each 8 x 8
  // block.) in a compressed form.
  memcpy(new_bits + BLAKE3_OUT_LEN, compressed_data, current_len);
  
  // Keep all the rest of the bits the same.
  memcpy(new_bits + BLAKE3_OUT_LEN + current_len, original_bits + BLAKE3_OUT_LEN + current_len,
	 original_bits_size - BLAKE3_OUT_LEN - current_len);
  
  for (i = y1; i < y2; i += 8)
    {
      /* Get row i through i+7 */
      gimp_pixel_rgn_get_row (&rgn_in,
                              row_arr[0],
                              x1, i,
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
                              row_arr[1],
                              x1, MIN (y2 - 1, i + 1),
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
			      row_arr[2],
                              x1, MIN (y2 - 1, i + 2),
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
			      row_arr[3],
                              x1, MIN (y2 - 1, i + 3),
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
			      row_arr[4],
                              x1, MIN (y2 - 1, i + 4),
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
			      row_arr[5],
                              x1, MIN (y2 - 1, i + 5),
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
			      row_arr[6],
                              x1, MIN (y2 - 1, i + 6),
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
                              row_arr[7],
                              x1, MIN (y2 - 1, i + 7),
                              x2 - x1);

      
      // Break up into 8x8 subblocks of pixels
      for (gint col_offset = 0; col_offset < channels * (x2 - x1); col_offset += 8){
        int x_block = col_offset / 8; // the column index of each block.
	int y_block = (i - y1) / 8; // the row index of each block
	int block_index = x_block + y_block * channels * width / 8; // there are channel * width / 8 blocks per row.

	// 2 bits per block are used. Gives the byte for the block.
	int original_bit_index = block_index / 4;

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
	// Get the byte with new_bits[original_index].
	//
	// Shift the byte so that the two bits we want are all the way
	// to the right with ">> (sub_block_index * 2))" (we multiply
	// sub_block_index by 2 because there are 2 bits per block.
	//
	// We use "& 3" to select the lowest two bits because other
	// bits are for other blocks.
	int new_value = (new_bits[original_bit_index] >> (sub_block_index * 2)) & 3;

	// Select the new values of the bits for each pixel.
	//
	// new_value has two bits. "&" it with 1 to get the lowest
	// bit. divide by 2 to get the next bit.
	// The order is arbitrary.
	int bit_1_p_alpha = (new_value & 1);
	int bit_alpha = (new_value / 2);

	// The 1_p_alpha pixel with the lowest bit set to 0.
	int row_1_p_alpha = (row_arr[y_bit_1_p_alpha_index][x_bit_1_p_alpha_index + col_offset] & 0xfe);

	// The alpha pixel with the lowest bit set to 0.
	int row_alpha = (row_arr[y_bit_alpha_index][x_bit_alpha_index + col_offset] & 0xfe);

	// Add the new bit value to the 1_p_alpha pixel.
	row_arr[y_bit_1_p_alpha_index][x_bit_1_p_alpha_index + col_offset] =
	  (row_1_p_alpha + bit_1_p_alpha);

	// Add the new bit value to the alpha pixel.
	row_arr[y_bit_alpha_index][x_bit_alpha_index + col_offset] =
	  (row_alpha + bit_alpha);
      }

      for (k = 0; k < 8; ++k){
	gimp_pixel_rgn_set_row (&rgn_out,
				row_arr[k],
				x1, i + k,
				x2 - x1);
      }
      
      if (i % 10 == 0)
	gimp_progress_update ((gdouble) (i - y1) / (gdouble) (y2 - y1));
    }

  for (i = 0; i < 8; ++i){
    g_free(row_arr[i]);
  }

  gimp_drawable_flush (drawable);
  gimp_drawable_merge_shadow (drawable->drawable_id, TRUE);
  gimp_drawable_update (drawable->drawable_id,
			x1, y1,
			x2 - x1, y2 - y1);
}
