#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>

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
  // guchar      *row1, *row2, *row3;
  guchar      *outrow;
  gint         width, height;

  gint u, v, x, y;

  guchar* row_arr[8];
  double* G_matrix_val[8]; // The matrix fot the DCT (Discrete Cosine Transform)
  double* G_prime_matrix_val[8];
  double* G_matrix_inverse_val[8];
  guchar* Q_matrix_val[8]; // The quantization matrix.
  // guchar* Q_inverse_matrix_val[8]; // The inverse of the quantization matrix
  int* B_matrix_val[8]; // The quantized DCT coefficient matrix

  int offset = 128;
  
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

  /* Initialise enough memory for row1, row2, row3, outrow */
  /* row1 = g_new (guchar, channels * (x2 - x1)); */
  /* row2 = g_new (guchar, channels * (x2 - x1)); */
  /* row3 = g_new (guchar, channels * (x2 - x1)); */
  /* row4 = g_new (guchar, channels * (x2 - x1)); */
  /* row5 = g_new (guchar, channels * (x2 - x1)); */
  /* row6 = g_new (guchar, channels * (x2 - x1)); */
  /* row7 = g_new (guchar, channels * (x2 - x1)); */
  /* row8 = g_new (guchar, channels * (x2 - x1)); */
  for (int i = 0; i < 8; ++i){
    row_arr[i] = g_new(guchar, channels * (x2 - x1));
    G_matrix_val[i] = g_new(double, channels * (x2 - x1));
    G_prime_matrix_val[i] = g_new(double, channels * (x2 - x1));
    G_matrix_inverse_val[i] = g_new(double, channels * (x2 - x1));
    Q_matrix_val[i] = g_new(guchar, channels * (x2 - x1));
    B_matrix_val[i] = g_new(int, channels * (x2 - x1));
  }

  printf("Crash point after initializing the matrices.\n");

  outrow = g_new (guchar, channels * (x2 - x1));

  // Create a quantization matrix

  Q_matrix_val[0][0] = 16;
  Q_matrix_val[0][1] = 11;
  Q_matrix_val[0][2] = 10;
  Q_matrix_val[0][3] = 16;
  Q_matrix_val[0][4] = 24;
  Q_matrix_val[0][5] = 40;
  Q_matrix_val[0][6] = 51;
  Q_matrix_val[0][7] = 61;
  Q_matrix_val[1][0] = 12;
  Q_matrix_val[1][1] = 12;
  Q_matrix_val[1][2] = 14;
  Q_matrix_val[1][3] = 19;
  Q_matrix_val[1][4] = 26;
  Q_matrix_val[1][5] = 58;
  Q_matrix_val[1][6] = 60;
  Q_matrix_val[1][7] = 55;
  Q_matrix_val[2][0] = 14;
  Q_matrix_val[2][1] = 13;
  Q_matrix_val[2][2] = 16;
  Q_matrix_val[2][3] = 24;
  Q_matrix_val[2][4] = 40;
  Q_matrix_val[2][5] = 57;
  Q_matrix_val[2][6] = 69;
  Q_matrix_val[2][7] = 56;
  Q_matrix_val[3][0] = 14;
  Q_matrix_val[3][1] = 17;
  Q_matrix_val[3][2] = 22;
  Q_matrix_val[3][3] = 29;
  Q_matrix_val[3][4] = 51;
  Q_matrix_val[3][5] = 87;
  Q_matrix_val[3][6] = 80;
  Q_matrix_val[3][7] = 62;
  Q_matrix_val[4][0] = 18;
  Q_matrix_val[4][1] = 22;
  Q_matrix_val[4][2] = 37;
  Q_matrix_val[4][3] = 56;
  Q_matrix_val[4][4] = 68;
  Q_matrix_val[4][5] = 109;
  Q_matrix_val[4][6] = 103;
  Q_matrix_val[4][7] = 77;
  Q_matrix_val[5][0] = 24;
  Q_matrix_val[5][1] = 25;
  Q_matrix_val[5][2] = 55;
  Q_matrix_val[5][3] = 64;
  Q_matrix_val[5][4] = 81;
  Q_matrix_val[5][5] = 104;
  Q_matrix_val[5][6] = 113;
  Q_matrix_val[5][7] = 92;
  Q_matrix_val[6][0] = 49;
  Q_matrix_val[6][1] = 64;
  Q_matrix_val[6][2] = 78;
  Q_matrix_val[6][3] = 87;
  Q_matrix_val[6][4] = 103;
  Q_matrix_val[6][5] = 121;
  Q_matrix_val[6][6] = 120;
  Q_matrix_val[6][7] = 101;
  Q_matrix_val[7][0] = 72;
  Q_matrix_val[7][1] = 92;
  Q_matrix_val[7][2] = 95;
  Q_matrix_val[7][3] = 98;
  Q_matrix_val[7][4] = 112;
  Q_matrix_val[7][5] = 100;
  Q_matrix_val[7][6] = 103;
  Q_matrix_val[7][7] = 99;


  for (i = y1; i < y2; i += 8)
    {
      /* Get row i through i+7 */
      printf("Crash point when getting row i.\n");
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
      // G_matrix_val_coord_6_6 = (1/4) * 1 * 1 * (coord 6,6) * cos((2 * 6 + 1) * 6 * pi() / 16) * cos((2 * 6 + 1) * 6 * pi() / 16)

      

      //Break up into 8x8 subblocks of pixels

      
      //for (j = x1; j < x2; j++){
	//if (j % 8 == 0){
	  for (u = 0; u < 8; ++u){ // u is the coordinate for the row_arr
	    for (v = 0; v <  8; ++v){
	      G_matrix_val[u][v] = 0;
		
	      for (x = 0; x <= 7; ++x){
		for (y = 0; y <= 7; ++y){
		  G_matrix_val[u][v] += (1.0/4) * alpha(u) * alpha(v) * (row_arr[x][y] - offset)
		    * cos((2 * x + 1) * u * M_PI / 16) * cos((2 * y + 1) * v * M_PI / 16);
		  /* if (u == 0 && v == 0){ */
		  /*   printf("%d %d %g %d %g %g \n", x, y, G_matrix_val[u][v], row_arr[x][y], alpha(u), alpha(v)); */
		  /* } */
		}
	      }
	    }
	  }
	  //	}
  //}

      printf("row_arr:\n");
      for (j = 0; j < 8; ++j){
	for (k = 0; k < 8; ++k){
	  printf("%d\t", row_arr[j][k]);
	}

	printf("\n");
      }

      printf("G_matrix_val:\n");
      for (j = 0; j < 8; ++j){
	for (k = 0; k < 8; ++k){
	  printf("%g\t", G_matrix_val[j][k]);
	}

	printf("\n");
      }

      // Compute the quantized DCT coefficient matrix
      for (j = 0; i < 8; ++i){
	for (k = 0; j < 8; ++j){
	  B_matrix_val[j][k] = round(G_matrix_val[j][k]/Q_matrix_val[j][k]);
	}
      }

      printf("B_matrix_val:\n");
      for (j = 0; j < 8; ++j){
	for (k = 0; k < 8; ++k){
	  printf("%d\t", B_matrix_val[j][k]);
	}

	printf("\n");
      }
      
      // Apply the paper to encode the watermark.
      // Compute the hash, and then encode it into the watermark (modify the matrix B).
      // Can use a fixed number as a hash. e.g. 1, 2, 3, 4, 5, 6, 7, 8, 9, ...
      // Encode it into the 6,6 component of the DCT.
      // Multiply the modified B element by element with Q (not regular multiply), and put the result in G'.

      for (j = 0; j < 8; ++j){
	for (k = 0; k < 8; ++k){
	  G_prime_matrix_val[j][k] = B_matrix_val[j][k] * Q_matrix_val[j][k];
	}
      }

      printf("G_prime_matrix_val:\n");
      for (j = 0; j < 8; ++j){
	for (k = 0; k < 8; ++k){
	  printf("%g\t", G_prime_matrix_val[j][k]);
	}

	printf("\n");
      }
      
      // Reverse the DCT on the new G' to get new values. Then you're done. (Multiply G' by the inverse of the DCT).
      
      // Compute the inverse matrix

      /* for (u = 0; u < 8; u++){ */
      /* 	for (v = 0; v < 8; v++){ */
      /* 	  G_matrix_inverse_val[u][v] = 0; */

      /* 	  for (x = 0; x < 8; x++){ */
      /* 	    for (y = 0; y < 8; ++y){ */
      /* 	      /\* G_matrix_inverse_val[u][v] += (1.0/4) * alpha(u) * alpha(v) * G_prime_matrix_val[x][y] *\/ */
      /* 	      /\* 	* cos((M_PI / 8) * (u + (1.0/2)) * x) * cos((M_PI / 8) * (v + (1.0/2)) * y); *\/ */

      /* 	      G_matrix_inverse_val[u][v] += (1.0/4) * alpha(u) * alpha(v) * G_matrix_val[x][y] */
      /* 		* cos((M_PI / 8) * (u + (1.0/2)) * x) * cos((M_PI / 8) * (v + (1.0/2)) * y); */
      /* 	      /\* if (u == 0 && v == 0){ *\/ */
      /* 	      /\* 	printf("%d %d %g \n", x, y, G_matrix_inverse_val[u][v]); *\/ */
      /* 	      /\* } *\/ */
      /* 	    } */
      /* 	  } */
      /* 	} */
      /* } */

      for (x = 0; x < 8; x++){
	for (y = 0; y < 8; y++){
	  G_matrix_inverse_val[x][y] = 0;

	  for (u = 0; u < 8; u++){
	    for (v = 0; v < 8; ++v){
	      /* G_matrix_inverse_val[u][v] += (1.0/4) * alpha(u) * alpha(v) * G_prime_matrix_val[x][y] */
	      /* 	* cos((M_PI / 8) * (u + (1.0/2)) * x) * cos((M_PI / 8) * (v + (1.0/2)) * y); */

	      G_matrix_inverse_val[x][y] += (1.0/4) * alpha(u) * alpha(v) * G_matrix_val[u][v]
		* cos((M_PI / 8) * (x + (1.0/2)) * u) * cos((M_PI / 8) * (y + (1.0/2)) * v);
	      /* if (u == 0 && v == 0){ */
	      /* 	printf("%d %d %g \n", x, y, G_matrix_inverse_val[u][v]); */
	      /* } */
	    }
	  }
	}
      }

      gint hash = 1;
      //G_matrix_val[6][6] = hash;

      
      printf("G_matrix_inverse_val:\n");
      for (j = 0; j < 8; ++j){
	for (k = 0; k < 8; ++k){
	  printf("%g\t", G_matrix_inverse_val[j][k] + offset);
	}

	printf("\n");
      }
      
      printf("Crash point after G_matrix_val\n");

      // inverse(Q_matrix_val[j][k], Q_inverse_matrix_val[j][k]);

      // multiplyMatrix(G_prime_matrix_val, G_matrix_inverse_val);
      
      // DCT_Coeff_matrix
      
      /* For each layer, compute the average of the nine
       * pixels */
      /* for (k = 0; k < channels; k++) */
      /* 	{ */
      /* 	  // outrow[channels * (j - x1) + k] = row2[channels * (j - x1) + k]; */
              
      /* 	  if (pixels_to_change[channels * (j - x1) * (i - y1) + k] == 1) { */
      /* 	    outrow[channels * (j - x1) + k] = 0; */
      /* 	  } */
      /* 	} */
      
      gimp_pixel_rgn_set_row (&rgn_out,
			      outrow,
			      x1, i,
			      x2 - x1);

      if (i % 10 == 0)
	gimp_progress_update ((gdouble) (i - y1) / (gdouble) (y2 - y1));
    }

  /* g_free (row1); */
  /* g_free (row2); */
  /* g_free (row3); */
  for (i = 0; i < 8; ++i){
    g_free(row_arr[i]);
    g_free(G_matrix_val[i]);
    g_free(Q_matrix_val[i]);
    g_free(B_matrix_val[i]);
  }
  g_free (outrow);

  printf(outrow);

  gimp_drawable_flush (drawable);
  gimp_drawable_merge_shadow (drawable->drawable_id, TRUE);
  gimp_drawable_update (drawable->drawable_id,
			x1, y1,
			x2 - x1, y2 - y1);
}
