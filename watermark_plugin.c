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

static void blur (GimpDrawable *drawable);

static gboolean
blur_dialog (GimpDrawable *drawable);

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

struct MyBlurVals {
  gint radius;
  gboolean preview;
};

/* Set up default values for options */
static struct MyBlurVals bvals; /* radius */
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

	/* switch (run_mode) */
	/*   { */
	/*   case GIMP_RUN_INTERACTIVE: */
	/*     /\* Get options last values if needed *\/ */
	/*     printf("Crash point at bvals.\n"); */
	/*     gimp_get_data ("plug-in-myblur", &bvals); */

	/*     /\* Display the dialog *\/ */
	/*     printf("Crash point in run.\n"); */
	/*     if (! blur_dialog (drawable)){ */
	/*       return; */
	/*     } */
	/*     printf("Crash point after blur_dialog.\n"); */
	/*     break; */

	/*   case GIMP_RUN_NONINTERACTIVE: */
	/*     if (nparams != 4) */
	/*       status = GIMP_PDB_CALLING_ERROR; */
	/*     if (status == GIMP_PDB_SUCCESS) */
	/*       bvals.radius = param[3].data.d_int32; */
	/*     break; */

	/*   case GIMP_RUN_WITH_LAST_VALS: */
	/*     /\*  Get options last values if needed  *\/ */
	/*     gimp_get_data ("plug-in-myblur", &bvals); */
	/*     break; */

	/*   default: */
	/*     break; */
	/*   } */
	
	gimp_progress_init ("My Blur...");

	/* Let's time blur
	 *
	 *   GTimer timer = g_timer_new time ();
	 */

	printf("Crash point just before blur.\n");	

	blur (drawable);

	/*   g_print ("blur() took %g seconds.\n", g_timer_elapsed (timer));
	 *   g_timer_destroy (timer);
	 */

	gimp_displays_flush ();
	gimp_drawable_detach (drawable);

	
	/*  Finally, set options in the core  */
	if (run_mode == GIMP_RUN_INTERACTIVE)
	  gimp_set_data ("plug-in-myblur", &bvals, sizeof (struct MyBlurVals));

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
blur (GimpDrawable *drawable)
{
  gint         i, j, k, channels;
  gint         x1, y1, x2, y2;
  GimpPixelRgn rgn_in, rgn_out;
  guchar      *row1, *row2, *row3;
  guchar      *outrow;
  gint         width, height;

  gimp_drawable_mask_bounds (drawable->drawable_id,
                             &x1, &y1,
                             &x2, &y2);
    width = x2 - x1;
    height = y2 - y1;
    
  channels = gimp_drawable_bpp (drawable->drawable_id);

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
  row1 = g_new (guchar, channels * (x2 - x1));
  row2 = g_new (guchar, channels * (x2 - x1));
  row3 = g_new (guchar, channels * (x2 - x1));
  outrow = g_new (guchar, channels * (x2 - x1));

  for (i = y1; i < y2; i++)
    {
      /* Get row i-1, i, i+1 */
      gimp_pixel_rgn_get_row (&rgn_in,
                              row1,
                              x1, MAX (y1, i - 1),
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
                              row2,
                              x1, i,
                              x2 - x1);
      gimp_pixel_rgn_get_row (&rgn_in,
                              row3,
                              x1, MIN (y2 - 1, i + 1),
                              x2 - x1);

      for (j = x1; j < x2; j++)
        {
          /* For each layer, compute the average of the nine
           * pixels */
          for (k = 0; k < channels; k++)
            {
              /* int sum = 0; */
              /* sum = row1[channels * MAX ((j - 1 - x1), 0) + k]           + */
              /*       row1[channels * (j - x1) + k]                        + */
              /*       row1[channels * MIN ((j + 1 - x1), x2 - x1 - 1) + k] + */
              /*       row2[channels * MAX ((j - 1 - x1), 0) + k]           + */
              /*       row2[channels * (j - x1) + k]                        + */
              /*       row2[channels * MIN ((j + 1 - x1), x2 - x1 - 1) + k] + */
              /*       row3[channels * MAX ((j - 1 - x1), 0) + k]           + */
              /*       row3[channels * (j - x1) + k]                        + */
              /*       row3[channels * MIN ((j + 1 - x1), x2 - x1 - 1) + k];  */
              outrow[channels * (j - x1) + k] = row2[channels * (j - x1) + k];
              // outrow[0] = 0;
            }

       }
      outrow[0] = 0;
      outrow[1] = 0;
      outrow[2] = 0;

      // outrow[0] = 0;
      // outrow = row;
       gimp_pixel_rgn_set_row (&rgn_out,
                               outrow,
                               x1, i,
                               x2 - x1);

       if (i % 10 == 0)
            gimp_progress_update ((gdouble) (i - y1) / (gdouble) (y2 - y1));
  }

  g_free (row1);
  g_free (row2);
  g_free (row3);
  g_free (outrow);

  gimp_drawable_flush (drawable);
  gimp_drawable_merge_shadow (drawable->drawable_id, TRUE);
  gimp_drawable_update (drawable->drawable_id,
                        x1, y1,
                        x2 - x1, y2 - y1);
}
