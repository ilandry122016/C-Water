#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include <jbig85.h>

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
  // fwrite(start, 1, len, (FILE *) file);
  
  return;
}

void multiplyMatrix(guchar** m1, guchar** m2)
{
  double result[8][8];

  printf("Resultant Matrix is:\n");

  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      result[i][j] = 0;

      for (int k = 0; k < 8; k++) {
	result[i][j] += m1[i][k] * m2[k][j];
      }

      printf("%d\t", result[i][j]);
    }

    printf("\n");
  }
}


/* Set up default values for options */
static struct MyWatermarkVals bvals; /* radius */
// bvals.radius = 3;

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
  printf("Crash point just before gtk_container_add.\n");
  // gtk_container_add (GTK_CONTAINER (GTK_DIALOG (dialog)->vbox), main_vbox);
  printf("Crash point just before gtk_widget_show.\n");
  gtk_widget_show (main_vbox);
  printf("Crash point just before gimp_drawable_preview_new.\n");

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

  /* Let's time watermark
   *
   *   GTimer timer = g_timer_new time ();
   */

  printf("Crash point just before watermark.\n");

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

  /*   g_print ("watermark() took %g seconds.\n", g_timer_elapsed (timer));
   *   g_timer_destroy (timer);
   */

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
  original_bits = g_new(guchar, channels * (width / 8) * (height / 8 ) / 4); // We need 2 bits per 8x8 block.

  gint max_difference = 0;

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
	
	  //printf("G %d %d %d %g \n", i, col_offset, (int)(G_matrix_val[6][6]), G_matrix_val[6][6]);
	  

      /* printf("row_arr:\n"); */
      /* for (j = 0; j < 8; ++j){ */
      /* 	for (k = 0; k < 8; ++k){ */
      /* 	  printf("%d\t", row_arr[j][k]); */
      /* 	} */

      /* 	printf("\n"); */
      /* } */

      /* printf("G_matrix_val:\n"); */
      /* for (j = 0; j < 8; ++j){ */
      /* 	for (k = 0; k < 8; ++k){ */
      /* 	  printf("%g\t", G_matrix_val[j][k]); */
      /* 	} */

      /* 	printf("\n"); */
      /* } */
      
      // Apply the paper to encode the watermark.
      // Compute the hash, and then encode it into the watermark (modify the matrix B).
      // Can use a fixed number as a hash. e.g. 1, 2, 3, 4, 5, 6, 7, 8, 9, ...
      // Encode it into the 6,6 component of the DCT.
      // Multiply the modified B element by element with Q (not regular multiply), and put the result in G'.

      for (j = 0; j < 8; ++j){
	for (k = 0; k < 8; ++k){
	  /* G_prime_matrix_val[j][k] = G_matrix_val[j][k]; */
	  G_prime_matrix_val[j][k] = 0;
	}
      }

      gint hash = 8; // Needs to be big enough to make a visible change
      G_prime_matrix_val[encode_u][encode_v] = hash;
      int x_block = col_offset / 8;
      int y_block = (i - y1) / 8;
      int block_index = x_block + y_block * channels * width / 8;
      int original_bit_index = block_index / 4;
      int sub_block_index = block_index & 3;

      int DCT_value = (int)(G_matrix_val[encode_u][encode_v]) & 3;
      original_bits[original_bit_index] = original_bits[original_bit_index] | (DCT_value << (sub_block_index * 2));
      printf("DCT_value: %d %d %d %d %d %d %d \n", x_block, y_block, block_index,
	     original_bit_index, sub_block_index, DCT_value, (int)(original_bits[original_bit_index]));

      /* printf("G_prime_matrix_val:\n"); */
      /* for (j = 0; j < 8; ++j){ */
      /* 	for (k = 0; k < 8; ++k){ */
      /* 	  printf("%g\t", G_prime_matrix_val[j][k]); */
      /* 	} */

      /* 	printf("\n"); */
      /* } */
      
      // Reverse the DCT on the new G' to get new values. Then you're done. (Multiply G' by the inverse of the DCT).

      for (x = 0; x < 8; x++){
	for (y = 0; y < 8; y++){
	  G_matrix_inverse_val[x][y] = 0;

	  for (u = 0; u < 8; u++){
	    for (v = 0; v < 8; ++v){
	      G_matrix_inverse_val[x][y] += (1.0/4) * alpha(u) * alpha(v) * G_prime_matrix_val[u][v]
		* cos((M_PI / 8) * (x + (1.0/2)) * u) * cos((M_PI / 8) * (y + (1.0/2)) * v);
	    }
	  }
	}

	// printf("%d %g \n", x, cos((M_PI / 8) * (x + (1.0/2)) * 6));
      }

      
      
     /*  printf("G_matrix_inverse_val:\n"); */
     /*  for (j = 0; j < 8; ++j){ */
     /* 	for (k = 0; k < 8; ++k){ */
     /* 	  printf("%g\t", G_matrix_inverse_val[j][k] + offset); */
     /* 	} */

     /* 	printf("\n"); */
     /*  } */

      for (k = 0; k < 8; ++k){
	for (j = 0; j < 8; ++j){
	  guchar high_bits = row_arr[k][j + col_offset] & 252;
	  guchar low_bits = row_arr[k][j + col_offset] & 3;
	  low_bits += G_matrix_inverse_val[j][k];
	  low_bits = low_bits % 4;
	  outrow[k][j + col_offset] = low_bits | high_bits;
	  // outrow[k][j + col_offset] = G_matrix_inverse_val[j][k] + row_arr[k][j + col_offset];
	  /* if (max_difference < abs((char)(outrow[k][j + col_offset] - row_arr[k][j + col_offset]))){ */
	  /*   max_difference = abs((char)(outrow[k][j + col_offset] - row_arr[k][j + col_offset])); */
	  /*   printf("%d %d %d %d %d %d \n", k, j, col_offset, max_difference, outrow[k][j + col_offset], row_arr[k][j + col_offset]); */
	  /* } */
	}
      }
      /* printf("row_arr[6][6] - outrow[6][6]: before: %d %d %d %g %g \n", col_offset, */
      /* 	     row_arr[6][6 + col_offset], outrow[6][6 + col_offset], */
      /* 	     G_matrix_inverse_val[6][6], */
      /* 	     row_arr[6][6 + col_offset] - (G_matrix_inverse_val[6][6] + offset)); */
      	

      }      
     /*  printf("Crash point after G_matrix_val\n"); */
      for (k = 0; k < 8; ++k){
	
	gimp_pixel_rgn_set_row (&rgn_out,
				outrow[k],
				x1, i + k,
				x2 - x1);
      }

      /* (1.0/4) * alpha(u) * alpha(v) * G_prime_matrix_val[u][v] */
      /* 		* cos((M_PI / 8) * (x + (1.0/2)) * u) * cos((M_PI / 8) * (y + (1.0/2)) * v); */
	    
      //printf("outrow[6][6]: after: %d \n", outrow[6][6]);

      if (i % 10 == 0)
	gimp_progress_update ((gdouble) (i - y1) / (gdouble) (y2 - y1));
    }

  printf("Crash point before freeing the memory.\n");
  for (i = 0; i < 8; ++i){
    printf("Crash point when freeing row_arr.\n");
    g_free(row_arr[i]);
    g_free(outrow[i]);
  }

  gimp_drawable_flush (drawable);
  gimp_drawable_merge_shadow (drawable->drawable_id, TRUE);
  gimp_drawable_update (drawable->drawable_id,
			x1, y1,
			x2 - x1, y2 - y1);

  unsigned char bitmap[15] = {
    /* 23 x 5 pixels, "JBIG" */
    0x7c, 0xe2, 0x38, 0x04, 0x92, 0x40, 0x04, 0xe2,
    0x5c, 0x44, 0x92, 0x44, 0x38, 0xe2, 0x38
  };
  struct jbg85_enc_state se;
  //int i;
  
  jbg85_enc_init(&se, 23, 5, output_bie, stdout);      /* initialize encoder */
  jbg85_enc_options(&se, JBG_TPBON, 0, -1);      /* clear JBG_VLENGTH option */
  for (i = 0; i < 5; i++) {
    /* encode line */
    jbg85_enc_lineout(&se, bitmap+i*3, bitmap+(i-1)*3, bitmap+(i-2)*3);
  }

  for (i = 0; i < 15; ++i){
    printf("%d \n",(int)(bitmap[i]));
  }
}
