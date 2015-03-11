#include "data_tools.h"


/* Initialize data structure */
int data_init (
        params_struct *params,
        data_struct *data
        ) {

  // Iterators
  int i_block=0, i=0;
  
  // Counter
  data->n_gals=0;

  // Read data from the params structure
  strcpy(data->out_filename, params->results_filename);

  // Allocate memory for the results arrays
  data->n_pairs = (long long*)calloc(params->n_separation_bins,sizeof(long long));

  // =========================
  // = Initialise the blocks =
  // =========================
  data->n_blocks = params->n_input_files;  
  data->blocks = (data_block_struct*)malloc(sizeof(data_block_struct)*data->n_blocks);

  // Print information on the statistic that we are going to compute
  if (params->modality == VINFALL_ISOLATION)
    printf("*** We shall now compute the INFALL velocity statistic with an ISOLATION criterion of r = %g, n_neighbours_max = %d\n",
      params->r_isolation,
      params->n_neighbours_max);
  else
    printf("*** We shall now compute the INFALL velocity statistic WITHOUT isolation criterion\n");

  
  // Initialise the blocks.  Most of the computations will be made in this loop.
  for(i=0; i<data->n_blocks; ++i) {

    if (data_block_init(params, &(data->blocks[i]), i) == _FAILURE_)
      return _FAILURE_;
      
  }

  // Compute the total number of galaxies from the dataset that we considered
  // in our computations.
  // N.B: the data->blocks[i].n_gals are filled only after data_block_init!
  for(i=0; i<data->n_blocks; ++i)
    data->n_gals += data->blocks[i].n_gals;

  
  // =========================================
  // = Calculate infall velocities statistic =
  // =========================================
  zelda_call(data_vinfall_stat (params, data));

  // Store results
  zelda_call(data_write_file(params, data));
    
  return _SUCCESS_;

}



/* Compute the infall velocities statistics related to full dataset,
 * starting from those already calculated for the single blocks
 */
int data_vinfall_stat (
        params_struct *params,
        data_struct *data
        ) {


  int i_bin=0, i_block=0;

  // Allocate memory for the result arrays (all with 'calloc' since they are
  // accumulators).  Their meaning is explained below in the double loop on
  // blocks & bins.

  // Pair-averaged quantities
  data->sum_vinfall = (double*)calloc(params->n_separation_bins,sizeof(double));
  data->avg_vinfall = (double*)calloc(params->n_separation_bins,sizeof(double));
  data->sum_vinfall_sq = (double*)calloc(params->n_separation_bins,sizeof(double));  
  data->variance_vinfall = (double*)calloc(params->n_separation_bins,sizeof(double));
  
  // Block-averaged quantities
  data->sum_vinfall_means = (double*)calloc(params->n_separation_bins,sizeof(double));
  data->avg_vinfall_means = (double*)calloc(params->n_separation_bins,sizeof(double));  
  data->sum_vinfall_means_sq = (double*)calloc(params->n_separation_bins,sizeof(double));
  data->variance_vinfall_means = (double*)calloc(params->n_separation_bins,sizeof(double));
  data->averaged_variance_vinfall = (double*)calloc(params->n_separation_bins,sizeof(double));


  // Accumulate the quantities needed to compute the statistics
  for(i_block=0; i_block<data->n_blocks; ++i_block) {

    for(i_bin=0; i_bin<params->n_separation_bins; ++i_bin) {

      // Add the number of pairs of the considered block for the considered bin
      // to the total of that bin
      data->n_pairs[i_bin] += data->blocks[i_block].n_pairs[i_bin];
      
      // Here we compute the total statistics relative to the full dataset, i.e.
      // to the ensemble of N galaxies (N = sum(n_i) over the blocks).
      // We still need to divide these quantities by the total number of pairs in
      // each bin, that is data->n_pairs[i_bin].  We will do this only outside
      // this loop, below.
      data->sum_vinfall[i_bin] += data->blocks[i_block].sum_vinfall[i_bin];
      data->sum_vinfall_sq[i_bin] += data->blocks[i_block].sum_vinfall_sq[i_bin];

      // Here instead we compute the block-averaged statistics, i.e. we take
      // as our random variable the (average) v_infall of each block.
      // Hence, the expectation values will be calculated using a population
      // of data->n_blocks blocks, rather than the full population of galaxies.
      // In this approach, We treat the averages of the relative infall velocities obtained
      // for each block as if they were independent measurements of that quantity.
      data->sum_vinfall_means[i_bin] += data->blocks[i_block].avg_vinfall[i_bin];
     
      data->sum_vinfall_means_sq[i_bin] += (data->blocks[i_block].avg_vinfall[i_bin]
                                            *data->blocks[i_block].avg_vinfall[i_bin]);
      
      // There also is a third kind of statistic, which is the average of the
      // block variances, that is avg(varV).  We have to think how it is related
      // to the other two variances.
      data->averaged_variance_vinfall[i_bin] += data->blocks[i_block].variance_vinfall[i_bin]/((double)(data->n_blocks));
      
    } // i_bin
  } // i_block
  
  data->n_pairs_total = 0;  
  // Compute the actual average (see comments above)
  for(i_bin=0; i_bin<params->n_separation_bins; ++i_bin) {

    // Here the average is over the total number of galaxies
    data->avg_vinfall[i_bin] = data->sum_vinfall[i_bin]/((double)(data->n_pairs[i_bin]));             // avg(V)
    data->variance_vinfall[i_bin] = data->sum_vinfall_sq[i_bin]/((double)(data->n_pairs[i_bin]))      // var(V)
                                    - data->avg_vinfall[i_bin]*data->avg_vinfall[i_bin];

    // Here the average is over the total number of blocks
    data->avg_vinfall_means[i_bin] = data->sum_vinfall_means[i_bin]/((double)(data->n_blocks));         // avg(avgV)
    
    data->variance_vinfall_means[i_bin] = data->sum_vinfall_means_sq[i_bin]/((double)(data->n_blocks))  // avg(varV)
                                          - data->avg_vinfall_means[i_bin]*data->avg_vinfall_means[i_bin];

    data->n_pairs_total += data->n_pairs[i_bin];
                                          
  }
  
  // Print some information
  if (params->verbose > 0) {
    
    zelda_test (data->n_gals==0, "stopping to prevent division by zero");
    
    printf("Full statistic: found %lld pairs with separations between r_min=%g and r_max=%g, using %d blocks (%.4g percent of possible pairs).\n",
      data->n_pairs_total, params->r_min, params->r_max, data->n_blocks,
      100.*data->n_pairs_total/((data->n_gals-1.)*data->n_gals/2.));

  }
  
  return _SUCCESS_;
  
}


/* Write to file the results for the whole dataset 
 */
int data_write_file (
          params_struct *params,
          data_struct *data
          ) {
  

  int i_bin;

  FILE *out_file = fopen(data->out_filename,"w");
  if (out_file==NULL) {
    printf("ERROR, %s: Cannot open file %s for writing.\n", __func__, data->out_filename);
    return _FAILURE_;
  }

  if (params->verbose > 0)
    printf("Storing results for the dataset in file %s.\n", data->out_filename);

  char labels[][128] = {
    "edges",
    "midpts",
    "npairs",
    "avgV",
    "varV",
    "sigmaV",
    "avgavgV",
    "varavgV",
    "sigmaavgV",
    "avgvarV"
  };

  fprintf(out_file, "# %28s %30s %30s %30s %30s %30s %30s %30s %30s %30s\n",
    labels[0], labels[1], labels[2], labels[3], labels[4], labels[5], labels[6], labels[7], labels[8], labels[9]);

  for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin) {
    fprintf(out_file, "%+30f %+30f %30lld %+30f %+30f %+30f %+30f %+30f %+30f %+30f\n", 
      params->r_edges[i_bin],
      params->r_midpoints[i_bin],
      data->n_pairs[i_bin],  
      data->avg_vinfall[i_bin],
      data->variance_vinfall[i_bin],
      sqrt(data->variance_vinfall[i_bin]), 
      data->avg_vinfall_means[i_bin],
      data->variance_vinfall_means[i_bin],
      sqrt(data->variance_vinfall_means[i_bin]),
      data->averaged_variance_vinfall[i_bin]
      );
  }

  fclose(out_file);

  return _SUCCESS_;

}



          
/* Read in the pre-computed files containing the linear infall velocities
 * of each block
 */
int data_read_vinfall_stat_from_files (
          params_struct *params,
          data_struct *data
          ) {
  
             
             
  return _SUCCESS_;
  
}
          
          
/* Free the memory associated to the structure */
int data_free (
        data_struct *data
        ) {
          
  int i_block;

  // Free blocks
  for(i_block=0; i_block<data->n_blocks; ++i_block) {
    if( data_block_free(&(data->blocks[i_block])) == _FAILURE_ )
      return _FAILURE_;
  }

  free(data->n_pairs);

  // Pair-averaged quantities
  free(data->sum_vinfall);
  free(data->avg_vinfall);
  free(data->sum_vinfall_sq);
  free(data->variance_vinfall);
  
  // Block-averaged quantities
  free(data->sum_vinfall_means);
  free(data->avg_vinfall_means);
  free(data->sum_vinfall_means_sq);
  free(data->variance_vinfall_means);
  free(data->averaged_variance_vinfall);
  
  return _SUCCESS_;
  
}
        
          
          
          
          
          
          
  