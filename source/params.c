#include "params.h"

/* Initialise params structure.  Basically, fill the structure using
 * the arguments passed to the function, create the separation binning
 * arrays, and build the filenames of the various blocks
 */
int params_init (
          params_struct *params
          ) {

  // Counters
  int i=0, j=0, k=0;

  // Set the sorting flag
  params->sorting = _FALSE_;
  if (params->modality == VINFALL_ISOLATION)
    params->sorting = _TRUE_;
    
  // If n_gals_per_block is something different from zero, then read the files only partially
  if (params->n_gals_per_block!=0)
    printf("IMPORTANT: Considering only the first %d galaxies of each block.\n", params->n_gals_per_block);
    
  // Number of files that will be read
  params->n_input_files = params->n_box_bins*params->n_box_bins*params->n_box_bins;


  // A few checks
  zelda_test(
    (params->binning_mode == LOG_BINNING) && ((params->r_min==0.) || (params->r_max==0.)),
    "A logarithmically spaced array cannot start or end with a zero!"
    );

  zelda_test(
    params->r_min > params->r_max,
    "Please provide r_min<r_max for the separation grid."
    );

  zelda_test(
    params->r_max > params->r_isolation,
    "Please provide r_max < r_isolation.  It does not make sense to ask for isolated pairs and then probe separations larger than the isolation criterion."
    );

  // Throw an error if 'n_neigbours_max==0'.  We cannot select any pair, if we do not allow galaxies
  // to have at least one other neighbouring galaxy.
  zelda_test(
    (params->n_neighbours_max==0) && (params->modality == VINFALL_ISOLATION),
    "you are after a 2-point statistic, please specify n_neigbours_max>0."
    );  


  // ==================
  // = Binning arrays =
  // ==================
  int i_bin=0;
  
  if (params->binning_mode == LOG_BINNING) {
    // Edges
    params->r_edges = log_space(params->r_min, params->r_max, params->n_separation_bins+1);
    // Logarithmic bin width
    params->r_bin_width = log(params->r_max/params->r_min)/((double)(params->n_separation_bins));
    zelda_test(params->r_bin_width==0, "stopping to avoid division by zero.");    
    // Midpoints
    zelda_alloc(params->r_midpoints, params->n_separation_bins*sizeof(double));
    for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin)
      params->r_midpoints[i_bin] = exp( log(params->r_edges[i_bin]) + 0.5*params->r_bin_width);
  }
  else if (params->binning_mode == LIN_BINNING) {
    // Edges
    params->r_edges = lin_space(params->r_min, params->r_max, params->n_separation_bins+1);
    // Bin width
    params->r_bin_width = (params->r_max-params->r_min)/((double)(params->n_separation_bins));
    zelda_test(params->r_bin_width==0, "stopping to avoid division by zero.");
    // Midpoints
    zelda_alloc(params->r_midpoints, params->n_separation_bins*sizeof(double));    
    for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin)
      params->r_midpoints[i_bin] = 0.5*(params->r_edges[i_bin+1] + params->r_edges[i_bin]);
  }
  else if (params->binning_mode == CUSTOM_BINNING) {
    // The bins edges and midpoints should have already been allocated and read from input.c.
    // The bin width is not needed.
  }
  

  
  // ==================================
  // = Build input & output filenames =
  // ==================================
  params->input_filenames = (char**)malloc(sizeof(char*)*params->n_input_files);
	if (params->store_blocks_results == _TRUE_)
  	params->output_filenames = (char**)malloc(sizeof(char*)*params->n_input_files);

  for(i=0; i<params->n_input_files; ++i) {
    params->input_filenames[i] = (char*)malloc(sizeof(char)*NAME_MAX);
		if (params->store_blocks_results == _TRUE_)
    	params->output_filenames[i] = (char*)malloc(sizeof(char)*NAME_MAX);    
  }
    
  // Create the output directory
	if (params->store_blocks_results == _TRUE_)
  	mkdir (dirname(params->output_format), 0777);
    
  for(i = 0; i < params->n_box_bins; ++i) {
    for(j = 0; j < params->n_box_bins; ++j) {
      for(k = 0; k < params->n_box_bins; ++k) {
      
        // The format of the input and output files is given by command line.
        // A typical choice would be respectively directory/gal_block_%d_%d_%d.dat
        // and directory/result_block_%d_%d_%d.dat
        int n_box_bins = params->n_box_bins;        
        int file_index = n_box_bins*n_box_bins*(i) + n_box_bins*j + k;
        sprintf(params->input_filenames[file_index], params->input_format, i+1, j+1, k+1);
				if (params->store_blocks_results == _TRUE_)
        	sprintf(params->output_filenames[file_index], params->output_format, i+1, j+1, k+1);

      }
    }
  }
  
  // =====================
  // = Print information =
  // =====================

  // Print the separation grid
  if (params->verbose > 1) {
    printf("We shall use the following grid in separation (Mpc/h, midpoints):\n");
    printf("*  ");
    for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin)
      printf("%.5g  ", params->r_midpoints[i_bin]);
    printf("\n");
  }

  // Print the data files that will be read
  if (params->verbose > 1)
    printf("We shall read in %d files/blocks.\n", params->n_input_files);

  if (params->verbose > 4) {
    for(i = 0; i < params->n_input_files; ++i)
      printf("%s\n", params->input_filenames[i]);
  }
  
  return _SUCCESS_;

}



/* Free the structure */
int params_free(
        params_struct *params
        ) {
          
  int i=0;

  free(params->r_edges);
  free(params->r_midpoints);
  
  for(i=0; i<params->n_input_files; ++i) {
    free(params->input_filenames[i]);
		if (params->store_blocks_results == _TRUE_)
    	free(params->output_filenames[i]);
  }
  
  return _SUCCESS_;
  
}




