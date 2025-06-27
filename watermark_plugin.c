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

/* // Function to get cofactor of mat[p][q] in cof[][]. n is */
/* // current dimension of mat[][] */
/* void getCof(int mat[N][N], int cof[N][N], int p, int q, int n) { */
/*     int i = 0, j = 0; */
/*     for (int row = 0; row < n; row++) { */
/*         for (int col = 0; col < n; col++) { */
/*             if (row != p && col != q) { */
/*                 cof[i][j++] = mat[row][col]; */
/*                 if (j == n - 1) { */
/*                     j = 0; */
/*                     i++; */
/*                 } */
/*             } */
/*         } */
/*     } */
/* } */

/* // Recursive function for finding determinant of matrix mat of dimension n */
/* int getDet(int mat[N][N], int n) { */
/*     if (n == 1) return mat[0][0]; */
/*     int det = 0; */
    
/*     int cof[N][N]; */
/*     int sign = 1; */
/*     for (int f = 0; f < n; f++) { */
/*         getCof(mat, cof, 0, f, n); */
/*         det += sign * mat[0][f] * getDet(cof, n - 1); */
/*         sign = -sign; */
/*     } */
/*     return det; */
/* } */

/* // Function to get adjoint of mat in adj */
/* void adjoint(int mat[N][N], double adj[N][N]) { */
/*     if (N == 1) { */
/*         adj[0][0] = 1; */
/*         return; */
/*     } */
    
/*     int sign = 1; */
/*     int cof[N][N]; */
/*     for (int i = 0; i < N; i++) { */
/*         for (int j = 0; j < N; j++) { */
/*             getCof(mat, cof, i, j, N); */
/*             sign = ((i + j) % 2 == 0) ? 1 : -1; */
/*             adj[j][i] = sign * getDet(cof, N - 1); */
/*         } */
/*     } */
/* } */

/* // Function to calculate and store inverse, returns 0 if matrix is singular */
/* int inverse(int mat[N][N], double inv[N][N]) { */
/*     int det = getDet(mat, N); */
/*     if (det == 0) { */
/*         printf("Singular matrix, can't find its inverse\n"); */
/*         return 0; */
/*     } */

/*     double adj[N][N]; */
/*     adjoint(mat, adj); */

/*     for (int i = 0; i < N; i++) */
/*         for (int j = 0; j < N; j++) */
/*             inv[i][j] = adj[i][j] / det; */

/*     return 1; */
/* } */

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
  
  /* for (gint i = y1; i < y2; i++) */
  /*   { */
  /*     for (gint j = x1; j < x2; j++) */
  /* 	{ */
  /* 	  // Get 8x8 subblocks of pixels */
  /* 	  // Set the value to 1 to change the pixel. */
  /* 	  // Set the value to 0 to leave the pixel alone. */
  /* 	  for (gint k = 0; k < channels; k++) { */
  /* 	    if (i % 8 == 0 & j % 8 == 0){ */
  /* 	      for (gint i_subblock = i; i_subblock < i_subblock + 8; i_subblock++){ */
  /* 		for (gint j_subblock = j; j_subblock < j_subblock + 8; j_subblock++){ */
  /* 		  // subblock_pixels[channels * (j - x1) * (i - y1) + k] = val_of_bit_image */
  /* 		} */
  /* 	      } */
  /* 	    } */
  /* 	    if (i == 6 && j == 6){ */

  /* 	      // G_matrix_val_coord_6_6 = (1/4) * 1 * 1 * (coord 6,6) * cos((2 * 6 + 1) * 6 * pi() / 16) * cos((2 * 6 + 1) * 6 * pi() / 16) */
  /* 	    } */
  /* 	    /\* if (i % 2 == 1 && j % 2 == 1){      *\/ */
  /* 	    /\*   pixels_to_change[channels * (j - x1) * (i - y1) + k] = 1; *\/ */
  /* 	    /\* } else { *\/ */
  /* 	    /\*   pixels_to_change[channels * (j - x1) * (i - y1) + k] = 0; *\/ */
  /* 	    /\* } *\/ */
  /* 	  } */
		
  /* 	} */
  /*   } */
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
  double* G_matrix_inverse_val[8];
  guchar* Q_matrix_val[8]; // The quantization matrix.
  // guchar* Q_inverse_matrix_val[8]; // The inverse of the quantization matrix
  guchar* B_matrix_val[8]; // The quantized DCT coefficient matrix
  
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
    G_matrix_val[i] = g_new(guchar, channels * (x2 - x1));
    G_matrix_inverse_val[i] = g_new(guchar, channels * (x2 - x1));
    Q_matrix_val[i] = g_new(guchar, channels * (x2 - x1));
    B_matrix_val[i] = g_new(guchar, channels * (x2 - x1));
  }

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
      
      for (j = x1; j < x2; j++){
	if (j % 8 == 0){
	  for (u = 0; u < 8; ++u){ // u is the coordinate for the row_arr
	    for (v = 0; v <  8; ++v){
	      G_matrix_val[u][v] = 0;
		
	      for (x = 0; x <= 7; ++x){
		for (y = 0; y <= 7; ++y){
		  G_matrix_val[u][v] = (1/4) * alpha(u) * alpha(v) * row_arr[x][y] * cos((2 * x + 1) * u * M_PI / 16) * cos((2 * y + 1) * v * M_PI / 16);
		}
	      }
	    }
	  }
	}
      }

      
      // Compute the quantized DCT coefficient matrix
      for (j = 0; i < 8; ++i){
	for (k = 0; j < 8; ++j){
	  B_matrix_val[j][k] = round(G_matrix_val[j][k]/Q_matrix_val[j][k]);
	}
      }

      // Apply the paper to encode the watermark.
      // Compute the hash, and then encode it into the watermark.
      // Can use a fixed number as a hash. e.g. 1, 2, 3, 4, 5, 6, 7, 8, 9, ...
      // Encode it into the 6,6 component of the DCT.
      // Reverse the DCT to get new values. Then you're done. (Multiply the DCT by it's inverse).
      
      // Compute the inverse matrix

      for (u = 0; u < 8; u++){
	for (v = 0; v < 8; v++){
	  G_matrix_inverse_val[u][v] = (1.0/4) * alpha(0) * alpha(0) * row_arr[0][0];

	  for (x = 1; x < 8; x++){
	    for (y = 1; y < 8; ++y){
	      G_matrix_inverse_val[u][v] += (1.0/4) * alpha(u) * alpha(v) * row_arr[x][y] * cos((M_PI / 8) * (u + (1/2)) * x) * cos((M_PI / 8) * (v + (1/2)) * y);
	    }
	  }
	}
      }

      gint hash = 1;
      G_matrix_val[6][6] = hash;

      // inverse(Q_matrix_val[j][k], Q_inverse_matrix_val[j][k]);

      multiplyMatrix(G_matrix_val, G_matrix_inverse_val);
      
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
