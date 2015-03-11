#ifndef _PARAMS_
#define _PARAMS_

#include "common.h"
#include "tools.h"

/* Params structure.  It contains all the general information that may be needed by all
 * other structures
 */
typedef struct {

  // File naming related
  FileName root;
  FileName input_format;
  FileName results_filename;
  int n_box_bins;
  int n_input_files;
  char ** input_filenames;
  
  // Cuts related
  int n_cuts;
  int * columns_to_cut;
  double * cuts;

  // Separation bins
  double r_min;
  double r_max;
  int n_separation_bins;
  
  // Enum structure that controls the desired type of binning
  enum BINNING{
        LIN_BINNING,
        LOG_BINNING,
        CUSTOM_BINNING
  } binning_mode;  
  
  // Edges & midpoints of the separation grid
  double *r_edges;
  double *r_midpoints;
  
  // The bin width is needed to automatically assign the pairs to the various bins.
  // It will be different whether the binning is linear or logarithmic.  For custom
  // binning, it is not needed (nor defined) at all.
  double r_bin_width;

  // Variable that defines how much of the simulation we want to process
  int n_gals_per_block;

  // Enum structure that controls which statistics are computed
  enum MODALITY{
        VINFALL,
        VINFALL_ISOLATION
  } modality;
    
  // Minimum isolation requested (used only if modality = VINFALL_ISOLATION),
  // in whatever distance units the input file is.
  double r_isolation;
  
  // Maximum number of neighbours allowed for a galaxy for it to be part
  // of a pair.  For the 'classical' isolation criterion used in Marinoni and
  // Buzzi (2010), this is 1.
  int n_neighbours_max;
  
  // Flag that controls whether the lines of each block should be sorted according
  // to their first column (usually the x position) before being processed.
  // This will be automatically enabled (==_TRUE_) when applying the isolation criterion
  // in calculating the infall velocity statistics.
  int sorting;

	// Should we store the intermediate results for each block?
	int store_blocks_results;
  FileName output_format;
  char ** output_filenames;    

  // Edge effects
  double side_of_box;
  int remove_galaxies_on_the_edges;

  // Verbosity flag  
  int verbose;
  
} params_struct;



/* Initialise params structure.  Basically, fill the structure using
 * the arguments passed to the function, create the separation binning
 * arrays, and build the filenames of the various blocks
 */
int params_init (
        params_struct *params
        );

/* Free the structure */
int params_free(
        params_struct *params
        );
        
#endif
