#include "zelda.h"

int main (int argc, char ** argv) {

  params_struct params;
  data_struct data;
  ErrorMsg error_message;

  // Cycle variable
  int i=0;

  // Read the input file
  zelda_call(input_init_from_arguments(argc, argv, &params, &data, error_message));

  printf("~~~~~ Value of my parameters: ~~~~~\n");

  printf("* General options\n");
  printf("\tparams.modality = %d\n", params.modality);
  
  printf("* Physical parameter\n");
  printf("\tparams.r_isolation = %g\n", params.r_isolation);
  printf("\tparams.n_neighbours_max = %d\n", params.n_neighbours_max);  
  
  printf("* File naming variables\n");
  printf("\tparams.root = %s\n", params.root);
  printf("\tparams.input_format = %s\n", params.input_format);
  printf("\tparams.results_filename = %s\n", params.results_filename);
	printf("\tparams.store_blocks_results = %d\n", params.store_blocks_results);
  printf("\tparams.output_format = %s\n", params.output_format);
  
  printf("* Cuts related parameters\n");
  printf("\tparams.n_cuts = %d\n", params.n_cuts);
  printf("\tparams.columns_to_cut = ");
  for(i=0; i<params.n_cuts; ++i)
    printf("%d  ", params.columns_to_cut[i]);
  printf("\n");
  printf("\tparams.cuts = ");
  for(i=0; i<2*params.n_cuts; ++i)
    printf("%.4g  ", params.cuts[i]);
  printf("\n");
  
  printf("* Removal of edges effects\n");
  printf("\tparams.side_of_box = %g\n", params.side_of_box);
  printf("\tparams.remove_galaxies_on_the_edges = %d\n", params.remove_galaxies_on_the_edges);  
  
  printf("* Binning\n");
  printf("\tparams.n_box_bins = %d\n", params.n_box_bins);
  printf("\tparams.r_min = %g\n", params.r_min);
  printf("\tparams.r_max = %g\n", params.r_max);
  printf("\tparams.n_separation_bins = %d\n", params.n_separation_bins);
  printf("\tparams.binning_mode = %d\n", params.binning_mode);
  if (params.binning_mode == CUSTOM_BINNING) {

    printf("\tparams.r_edges = ");
    for(i=0; i<params.n_separation_bins+1; ++i)
      printf("%.4g  ", params.r_edges[i]);
    printf("\n  ");

    printf("\tparams.r_midpoints = ");
    for(i=0; i<params.n_separation_bins; ++i)
      printf("%.4g  ", params.r_midpoints[i]);
    printf("\n");
  }

  
  printf("* Debug\n");
  printf("\tparams.n_gals_per_block = %d\n", params.n_gals_per_block);  
  printf("\tparams.verbose = %d\n", params.verbose);    
  
  return _SUCCESS_;
  
}