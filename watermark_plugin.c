#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
// #include <jbig85.h>
#include <jbig.h>
#include "blake3.h"
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

static void query (void);
static void run (const gchar *name,
		 gint nparams,
		 const GimpParam *param,
		 gint *nreturn_vals,
		 GimpParam **return_vals);

void gimp_pixel_rgn_get_pixel (GimpPixelRgn *pr,
			       guchar *buf,
			       gint x,
			       gint y);
void gimp_pixel_rgn_get_row (GimpPixelRgn *pr,
			     guchar *buf,
			     gint x,
			     gint y,
			     gint width);
void gimp_pixel_rgn_get_col (GimpPixelRgn *pr,
			     guchar *buf,
			     gint x,
			     gint y,
			     gint width);
void gimp_pixel_rgn_get_rect (GimpPixelRgn *pr,
			      guchar *buf,
			      gint x,
			      gint y,
			      gint width,
			      gint height);
void gimp_pixel_rgn_set_pixel (GimpPixelRgn *pr,
			       const guchar *buf,
			       gint x,
			       gint y);
void gimp_pixel_rgn_set_row (GimpPixelRgn *pr,
			     const guchar *buf,
			     gint x,
			     gint y,
			     gint width);
void gimp_pixel_rgn_set_col (GimpPixelRgn *pr,
			     const guchar *buf,
			     gint x,
			     gint y,
			     gint height);
void gimp_pixel_rgn_set_rect (GimpPixelRgn *pr,
			      const guchar *buf,
			      gint x,
			      gint y,
			      gint width,
			      gint height);

static void watermark (GimpDrawable *drawable, guchar *pixels_to_change, gint lower_limit_x, gint lower_limit_y, gint upper_limit_x, gint upper_limit_y, gint channels);

static gboolean
watermark_dialog (GimpDrawable *drawable);

static void init_mem (guchar ***row,
		      guchar **outrow,
		      gint num_bytes);
static void process_row (guchar **row,
			 guchar *outrow,
			 gint x1,
			 gint y1,
			 gint width,
			 gint height,
			 gint channels,
			 gint i);
static void shuffle (GimpPixelRgn *rgn_in,
		     guchar **row,
		     gint x1,
		     gint y1,
		     gint width,
		     gint height,
		     gint ypos);

// global variables
unsigned char compressed_data[1000];
size_t current_len = 0;

// helper function
static double alpha(gint i){
  if (i == 0){
    return 1.0/sqrt(2);
  }
  return 1;
}

struct MyWatermarkVals {
  gint radius;
  gboolean preview;
};

void output_bie(unsigned char *start, size_t len, void *file)
{
  for (int i = 0; i < len; ++i){
    compressed_data[current_len + i] = start[i];
  }
  current_len += len;
  
  return;
}

/* Set up default values for options */
static struct MyWatermarkVals bvals; /* radius */

/* The radius is still a constant, we'll change that when the
 * graphical interface will be built. */
static gint radius = 3;

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
	
  gimp_ui_init("mywatermark", FALSE);
  GtkWidget *main_vbox;
  main_vbox = gtk_vbox_new (FALSE, 6);
  gtk_widget_show (main_vbox);

  /* Setting mandatory output values */
  *nreturn_vals = 1;
  *return_vals  = values;
  bvals.radius = 3;

  values[0].type = GIMP_PDB_STATUS;
  values[0].data.d_status = status;

  /* Getting run_mode - we won't display a dialog if 
   * we are in NONINTERACTIVE mode */
  run_mode = param[0].data.d_int32;

  /*  Get the specified drawable  */
  drawable = gimp_drawable_get (param[2].data.d_drawable);
	
  gimp_progress_init ("My Watermark...");

  gint x1, y1, x2, y2;
	
  gimp_drawable_mask_bounds (drawable->drawable_id,
                             &x1, &y1,
                             &x2, &y2);

  guchar *pixels_to_change;
  gint channels = gimp_drawable_bpp (drawable->drawable_id);
  pixels_to_change = g_new(guchar, channels * (x2 - x1) * (y2 - y1));

  GimpPixelRgn rgn_in, rgn_out;
  guchar      *row1, *row2, *row3;
  guchar      *outrow;
  gint         width, height;

  guchar *subblock_pixels;
  
  watermark (drawable, pixels_to_change, x1, x2, y1, y2, channels);

  gimp_displays_flush ();
  gimp_drawable_detach (drawable);

	
  /*  Finally, set options in the core  */
  if (run_mode == GIMP_RUN_INTERACTIVE)
    gimp_set_data ("plug-in-mywatermark", &bvals, sizeof (struct MyWatermarkVals));

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
			  "plug-in-my-watermark",
			  "my-watermark",
			  "Displays \"My watermark\" in a dialog",
			  "Isaac Landry",
			  "Copyright Isaac Landry",
			  "2004",
			  "_My watermark...",
			  "RGB*, GRAY*",
			  GIMP_PLUGIN,
			  G_N_ELEMENTS (args), 0,
			  args, NULL);

  gimp_plugin_menu_register ("plug-in-my-watermark",
			     "<Image>/Filters/Misc");
}

static void
watermark(GimpDrawable *drawable, guchar *pixels_to_change, gint lower_limit_x, gint lower_limit_y, gint upper_limit_x, gint upper_limit_y, gint channels)
{
  gint         i, j, k;
  
  gint x1 = lower_limit_x;
  gint y1 = lower_limit_y;
  gint x2 = upper_limit_x;
  gint y2 = upper_limit_y;
  
  GimpPixelRgn rgn_in, rgn_out;
  guchar      *outrow[8];
  gint         width, height;

  gint u, v, x, y;

  int encode_u = 6;
  int encode_v = 6;

  guchar* row_arr[8];
  double G_matrix_val[8][8]; // The matrix fot the DCT (Discrete Cosine Transform)
  double G_prime_matrix_val[8][8];
  double G_matrix_inverse_val[8][8];
  int offset = 128; // To map the values from 0-255 to -128-127
  guchar* original_bits;
  
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
    outrow[i] = g_new(guchar, channels * (x2 - x1));
  }

  size_t original_bits_size = channels * (width / 8) * (height / 8 ) / 4;
  original_bits = g_new(guchar, original_bits_size); // We need 2 bits per 8x8 block. One is for the signal, and the other is for the modulator.

  gint max_difference = 0;

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
      }


      // Create A DCT matrix for the rows.
      
      

      //Break up into 8x8 subblocks of pixels

      
      for (gint col_offset = 0; col_offset < channels * (x2 - x1); col_offset += 8){
	for (u = 0; u < 8; ++u){ // u is the coordinate for the row_arr
	  for (v = 0; v <  8; ++v){
	    G_matrix_val[u][v] = 0;
		
	    for (x = 0; x <= 7; ++x){
	      for (y = 0; y <= 7; ++y){
		G_matrix_val[u][v] += (1.0/4) * alpha(u) * alpha(v) * (row_arr[y][col_offset + x] - offset)
		  * cos((2 * x + 1) * u * M_PI / 16) * cos((2 * y + 1) * v * M_PI / 16);
	      }
	    }
	  }
	}
			
	int x_block = col_offset / 8;
	int y_block = (i - y1) / 8;
	int block_index = x_block + y_block * channels * width / 8;
	int original_bit_index = block_index / 4;
	int sub_block_index = block_index & 3;

	// Get the lowest 2 integer bits of G[6][6].
	int DCT_value = (int)(G_matrix_val[encode_u][encode_v]) & 3;
	original_bits[original_bit_index] = original_bits[original_bit_index] | (DCT_value << (sub_block_index * 2));
	     
      }      

      if (i % 10 == 0)
	gimp_progress_update ((gdouble) (i - y1) / (gdouble) (y2 - y1));
    }

   // Finalize the hash. BLAKE3_OUT_LEN is the default output length, 32 bytes.
  uint8_t blake3_hash[BLAKE3_OUT_LEN];
  blake3_hasher_finalize(&hasher, blake3_hash, BLAKE3_OUT_LEN);

  gimp_drawable_flush (drawable);
  gimp_drawable_merge_shadow (drawable->drawable_id, TRUE);
  gimp_drawable_update (drawable->drawable_id,
			x1, y1,
			x2 - x1, y2 - y1);

  unsigned char *bitmaps[1] = {original_bits};
  struct jbg_enc_state se;
 
  jbg_enc_init(&se, (channels * width / 8) / 4, height / 8, 1, bitmaps, 
	       output_bie, stdout);              /* initialize encoder */
  jbg_enc_out(&se);                                    /* encode image */
  jbg_enc_free(&se);                    /* release allocated resources */

  guchar* new_bits = g_new(guchar, original_bits_size);

  memcpy(new_bits, blake3_hash, BLAKE3_OUT_LEN);
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
	// Create A DCT matrix for the 8x8 block.
	for (u = 0; u < 8; ++u){ // u is the coordinate for the row_arr
	  for (v = 0; v <  8; ++v){
	    G_matrix_val[u][v] = 0;
		
	    for (x = 0; x <= 7; ++x){
	      for (y = 0; y <= 7; ++y){
		G_matrix_val[u][v] += (1.0/4) * alpha(u) * alpha(v) * (row_arr[y][col_offset + x] - offset)
		  * cos((2 * x + 1) * u * M_PI / 16) * cos((2 * y + 1) * v * M_PI / 16);
	      }
	    }
	  }
	}
			
	int x_block = col_offset / 8;
	int y_block = (i - y1) / 8;
	int block_index = x_block + y_block * channels * width / 8;
	int original_bit_index = block_index / 4;
	int sub_block_index = block_index & 3;

	// We change values of the original image such that the DCT changes to what we want.
	// Get the lowest 2 integer bits of G[6][6].
	int DCT_value = (int)(G_matrix_val[encode_u][encode_v]) & 3;
	original_bits[original_bit_index] = original_bits[original_bit_index] | (DCT_value << (sub_block_index * 2));
	     
      }      

      if (i % 10 == 0)
	gimp_progress_update ((gdouble) (i - y1) / (gdouble) (y2 - y1));
    }

  for (i = 0; i < 8; ++i){
    g_free(row_arr[i]);
    g_free(outrow[i]);
  }
}
