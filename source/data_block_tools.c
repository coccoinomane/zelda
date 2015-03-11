#include "data_block_tools.h"


/* Initialise block structure */
int data_block_init (
    params_struct *params,
    data_block_struct *block,
    int block_id
    ) {

      
  // ==================
  // =   Read in file =
  // ==================

  // Assign id and name to the block according to the argument passed to the function
  block->id = block_id;
  strcpy(block->in_filename, params->input_filenames[block->id]);
	if (params->store_blocks_results == _TRUE_)
  	strcpy(block->out_filename, params->output_filenames[block->id]);
  
  // Find the position of the block in the 3D grid in x,y,z.
  // TODO: use a better way using n % n_bins etc.
  if (params->n_input_files==1)
    // If the data is not split, then we cannot infer i_x,i_y,i_z from the filename,
    // and we should set them by hand to 1
    block->i_x = block->i_y = block->i_z = 1;
  else
    sscanf(block->in_filename, params->input_format, &block->i_x, &block->i_y, &block->i_z);

  // Quantity needed to remove the edges effects
  block->side_length = params->side_of_box/(double)params->n_box_bins;
  
  // Print some information about the file that is being processed
  if (params->verbose > 1)
    printf("Processing block %d, corresponding to the file %s.\n", block->id+1, block->in_filename);

  // Read in the file associated to the block, thus filling the position and
  // velocity vectors and block->n_gals from the block file
  zelda_call(data_block_read_file(params, block));

     
  // ===================
  // = Allocate arrays =
  // ===================
  // Now that we know how many galaxies are in the block, we can allocate memory
  // for the remaining temporary & result arrays
  
  // Number of pairs that fall in a certain bin
  block->n_pairs = (long long*)calloc(params->n_separation_bins,sizeof(long long));

  // =======================================
  // = Calculate infall velocity statistic =
  // =======================================
  if (params->modality == VINFALL_ISOLATION) {
		zelda_call (data_block_vinfall_stat_with_isolation(params, block));
  } else {
		zelda_call (data_block_vinfall_stat(params, block));  
  }
  
  // Store results relative to the block, but only if requested
	if (params->store_blocks_results == _TRUE_)
	  if (data_block_write_file(params, block) == _FAILURE_)
	    return _FAILURE_;
  
  // ===============
  // = Free memory =
  // ===============
  data_block_free_workspace(block);

  return _SUCCESS_;
  
}





/* Calculate statistics considering all pairs of galaxies */
int data_block_vinfall_stat (
          params_struct *params,
          data_block_struct *block
          ) {


  long long i_gal=0, j_gal=0;
  int i_bin=0;

  // Allocate memory for the result arrays
  block->sum_vinfall = (double*)calloc(params->n_separation_bins,sizeof(double));
  block->avg_vinfall = (double*)calloc(params->n_separation_bins,sizeof(double));
  block->sum_vinfall_sq = (double*)calloc(params->n_separation_bins,sizeof(double));
  block->variance_vinfall = (double*)calloc(params->n_separation_bins,sizeof(double));

  // This variable may be useful if we are dealing with log spaced separation bins.
  // We compute it outside the main loop.
  double log_r_min;
  if (params->binning_mode == LOG_BINNING)
    log_r_min = log(params->r_min);

  // ==========================
  // = Main cycle on galaxies =
  // ==========================
  // Accumulator to count the pairs with zero distance
  long long n_same_position_gals=0;
  
  if (params->verbose > 2)
    printf("%s: Entering the main loop.\n", __func__);

  // BEGINNING OF PARALLEL LOOP
  # pragma omp parallel for private(i_gal,j_gal,i_bin) schedule(dynamic,1000)
      
  for (i_gal=0; i_gal < block->n_gals-1; ++i_gal) {
      
    for (j_gal=i_gal+1; j_gal<block->n_gals; ++j_gal) {

      double distance = get_distance(block->xx[i_gal],block->xx[j_gal],
                                     block->yy[i_gal],block->yy[j_gal],
                                     block->zz[i_gal],block->zz[j_gal]);

      // If a point-like pair is found (e.g. a pair made of just a galaxy that is repeated in the input
      // file), skip it and issue a warning message.
      // Note that this is not enough to prevent wierd results to show up.  For example, if you have two
      // repeated galaxies, you may have two pairs gal1-gal2 which are exactly the same.  This actually 
      // happens!  You can see it by looking for bins with only two pairs, but that have zero variance.
      // We check for this possibility just after this main cycle.
      if (distance<=0.) {
        #pragma omp atomic
        n_same_position_gals += 1;
        if (params->verbose > 4)
          WARNINGPRINT("Found a <= zero distance for the pair (%lld,%lld); distance = %lf\n", i_gal, j_gal, distance);
        continue;
      }

      // Check that the pair is in the requested range, otherwise skip it.
      if ((distance < params->r_edges[0]) || (distance >= params->r_edges[params->n_separation_bins]))
        continue;
        
      // Assign the pair to its separation bin.
      if (params->binning_mode == LIN_BINNING)
        i_bin = (int)((distance-params->r_min)/params->r_bin_width);
        
      else if (params->binning_mode == LOG_BINNING)
        i_bin = (int)((log(distance) - log_r_min)/params->r_bin_width);
        
      else if (params->binning_mode == CUSTOM_BINNING) {
        for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin)
          if ( (distance >= params->r_edges[i_bin]) && (distance < params->r_edges[i_bin+1]) )
            continue;
      }

      double relative_vinfall =  get_relative_infall_vel(block->xx[i_gal],block->xx[j_gal],
                                                        block->yy[i_gal],block->yy[j_gal],
                                                        block->zz[i_gal],block->zz[j_gal],
                                                        block->vxx[i_gal],block->vxx[j_gal],
                                                        block->vyy[i_gal],block->vyy[j_gal],
                                                        block->vzz[i_gal],block->vzz[j_gal],
                                                        distance);


      // Compute the quantities we need
      #pragma omp atomic
      block->n_pairs[i_bin] += 1;
      #pragma omp atomic
      block->sum_vinfall[i_bin] += relative_vinfall;
      #pragma omp atomic
      block->sum_vinfall_sq[i_bin] += (relative_vinfall*relative_vinfall);
        
      //// Debug.  Comment it if you need more speed.
      // if (params->verbose > 2) 
      //   printf("pair (%d,%d): distance = %5f, relative_vinfall = %5f\n", i_gal, j_gal, distance, relative_vinfall);
      
      //// Check our binning technique. Comment only if you need more speed.
      //// Double check on the positioning of pairs in the separation bins.
      //// When dealing with a logarithmically spaced grid, every now and then
      //// a pair fails this check, because of the < / <= and > / >= ambiguity.
      //// I guess it is an issue of precision.  We neglect this effect since
      //// it will "misplace" just 1 pair out of > 10 millions.
      // int i_temp=0;
      // for (i_temp=0; i_temp<params->n_separation_bins; ++i_temp) {
      //   // NB: In order to compare with the above bin assignation method, we need
      //   // to assign pairs lying on the division between two adjacent bins to the
      //   // LEFT bin.
      //   if ( (distance>=params->r_edges[i_temp]) && (distance<params->r_edges[i_temp+1]) &&
      //        (i_temp != i_bin) ) {
      //     
      //     WARNINGPRINT("Pair (%d,%d) (distance %lf, log(distance/r_min) %lf, params->r_bin_width %lf) of block %d_%d_%d was put in bin %d instead that in bin %d\n",
      //       i_gal, j_gal, distance, log(distance)-log_r_min, params->r_bin_width, block->i_x, block->i_y, block->i_z, i_bin, i_temp);
      // 
      //   }
      // }      
        
    } // End of loop on j_gal
                  
  } // End of loop on i_gal
  

  // Compute the statistics for the block
  for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin) {

    // Average
    block->avg_vinfall[i_bin] = block->sum_vinfall[i_bin]/((double)(block->n_pairs[i_bin]));

    // The variance is <V_j^2> - <V_j>^2 , where j is a pair's index
    block->variance_vinfall[i_bin] = block->sum_vinfall_sq[i_bin]/((double)(block->n_pairs[i_bin]))
                                     - block->avg_vinfall[i_bin]*block->avg_vinfall[i_bin];
   
    if((block->variance_vinfall[i_bin] == 0.) && (block->n_pairs[i_bin]>1) && params->verbose>0)
      WARNINGPRINT("The bin %d for the block %d has zero variance but %lld pairs!\n", i_bin, block->id, block->n_pairs[i_bin]);
    
  }
  
  // Issue a warning if some point-like pairs were found
  if (params->verbose>0 && n_same_position_gals!=0)
    WARNINGPRINT("%lld point-like pairs of galaxies %s found and excluded!\n",
      n_same_position_gals, n_same_position_gals>1 ? "were":"was");
    
  return _SUCCESS_;
  
}



/* Compute the pairwise velocity statistics for isolated galaxy pairs.
 * The degree of isolation is controlled by 'params->r_isolation' and
 * 'params->n_neighbours_max'.
 */
int data_block_vinfall_stat_with_isolation (
           params_struct *params,
           data_block_struct *block
           )
{

  // If we are looking for pairs that are completely isolated, then we can use
  // the more specific function 'data_block_vinfall_stat_with_total_isolation'.  
  // The execution speed is more or less the same, but it uses a different
  // algorithm.   Uncomment the following 'if' block to use the specific
  // function.
  // if (params->n_neighbours_max==1) {
  //   zelda_call (data_block_vinfall_stat_with_total_isolation(params, block));
  //   return _SUCCESS_;
  // }

  long long i_gal=0, j_gal=0;
  long long n_pairs=0, n_same_position_gals=0;
  int i_bin=0;

  // Determine the number of neighbours for each galaxy.  This is a crucial step to
  // apply the algorithm that will be employed in the main cycle, below.  It is likely
  // the step that takes most of the execution time.
  if (params->verbose > 2)
    printf("Computing the number of neighbours of each galaxy... \n");

  zelda_call(
    data_block_neighbours (params, block)
    );

  // Allocate memory for the result arrays
  block->sum_vinfall = (double*)calloc(params->n_separation_bins,sizeof(double));
  block->avg_vinfall = (double*)calloc(params->n_separation_bins,sizeof(double));
  block->sum_vinfall_sq = (double*)calloc(params->n_separation_bins,sizeof(double));
  block->variance_vinfall = (double*)calloc(params->n_separation_bins,sizeof(double));

  // This variable may be useful if we are dealing with log spaced separation bins.
  // We compute it outside the main loop.
  double log_r_min;
  if (params->binning_mode == LOG_BINNING)
    log_r_min = log(params->r_min);

  if (params->verbose > 3)
    printf("%s: Entering the main loop.\n", __func__);

  // ==========================
  // = Main cycle on galaxies =
  // ==========================
  //  Beginning of parallel region
  # pragma omp parallel for private(i_gal, j_gal, i_bin) schedule(dynamic,1000)

  // We shall exploit the fact that the galaxies have been sorted in the x-direction
  // to allow a faster computation.
  // For each galaxy i_gal, we shall only consider only
  // the other galaxies whose x-coordinate is close enough to the x-coordinate
  // of 'i_gal'.  Since the galaxies are already sorted according to x, this can
  // be done very quickly by splitting the loop on the j_gals
  // in a forward and a backward loop..
  for (i_gal=0; i_gal < block->n_gals-1; ++i_gal) {

    // If 'i_gal' has less than zero neighbours, there must be an error
    if (block->n_neighbours[i_gal]<0 && block->n_neighbours[i_gal]!=_BADFLAG_)
      ERRORPRINT("Galaxy %lld has a negative number (%lld) of neighbours!\n", i_gal, block->n_neighbours[i_gal]);

    // =======================================
    // = Cut in the # of neighbours of i_gal =
    // =======================================
    // Skip the galaxy if it has zero neighbours, or too many neighbours.
    if (block->n_neighbours[i_gal]==0 || block->n_neighbours[i_gal]==_BADFLAG_)
      continue;

    // We shall now start the loop on all other galaxies, i.e. on 'j_gal'.  We shall keep track
    // of the number of neighbours we find for 'i_gal', and stop the loop when this number is 
    // equal to the total number of pairs for 'i_gal'.
    long long n_igal_neighbours_so_far=0;

    for (j_gal=i_gal+1; j_gal<block->n_gals; ++j_gal) {

      // =========================
      // = Cut in the x-position =
      // =========================
      // Distance in x between the pair
      double delta_x = block->xx[j_gal]-block->xx[i_gal];

      // Check that x is sorted
      if (delta_x < 0)
        WARNINGPRINT("The input file %s has not been sorted as it should.\n",block->in_filename);

      // Stop looking for pairs when the delta_x is bigger than the isolation criterion.
      if( delta_x > params->r_isolation )
        break;


      // =======================================
      // = Cut in the # of neighbours of j_gal =
      // =======================================
      // Skip the pair if the second galaxy has too many neighbours.  This condition
      // is necessary because we want to get rid of any 'bad_flagged' galaxy (see
      // comment above).
      if(block->n_neighbours[j_gal]==_BADFLAG_)
        continue;


      // ===================
      // = Cut in distance =
      // ===================
      double distance = get_distance(block->xx[i_gal],block->xx[j_gal],
                                  block->yy[i_gal],block->yy[j_gal],
                                  block->zz[i_gal],block->zz[j_gal]);

      // Skip this pair if it is far away.
      if (distance > params->r_isolation)
        continue;


      // Perform a check on point-like pairs, i.e. repeated lines in the input files
      if(distance==0.) {
        if (params->verbose>3) {
          WARNINGPRINT("Pair (%lld,%lld) is point-like!\n",i_gal+1,j_gal+1);
          if( (block->vxx[i_gal]!=block->vxx[j_gal]) )
            printf("         (the velocity is different though!)\n");
        }
        #pragma omp atomic
        ++n_same_position_gals;
        continue;
      }

      // At this point we can be reasonably sure that (i_gal,j_gal) is a proper pair.
      #pragma omp atomic
      ++n_pairs;
      if (params->verbose > 3) {
        printf("Galaxies (%lld,%lld)-1 are matched (distance=%g).  Now we have %lld pairs.\n", i_gal+1, j_gal+1, distance, n_pairs);
        printf("    x=%20g, y=%20g, z=%20g\n", block->xx[i_gal], block->yy[i_gal], block->zz[i_gal]);
        printf("    x=%20g, y=%20g, z=%20g\n", block->xx[j_gal], block->yy[j_gal], block->zz[j_gal]);        
      }

      // Small consistency check
      if (block->n_neighbours[j_gal]==0)
        ERRORPRINT("How is it possible that 'j_gal' has 0 neighbours if we just paired it with 'i_gal'?\n","");
      

      // ===========================
      // = Compute the statistics =
      // ==========================

      // Check that the pair is in the requested range, otherwise skip it.
      if ((distance < params->r_edges[0]) || (distance >= params->r_edges[params->n_separation_bins]))
        continue;
        
      // Assign the pair to its separation bin.
      if (params->binning_mode == LIN_BINNING)
        i_bin = (int)((distance-params->r_min)/params->r_bin_width);

      else if (params->binning_mode == LOG_BINNING)
        i_bin = (int)((log(distance) - log_r_min)/params->r_bin_width);

      else if (params->binning_mode == CUSTOM_BINNING) {
       
        for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin)
          if ( (distance >= params->r_edges[i_bin]) && (distance < params->r_edges[i_bin+1]) )
            break;
      }

      if (params->verbose > 3)
        printf("Pair (%lld,%lld)-1 has been put into i_bin=%d.\n", i_gal+1, j_gal+1, i_bin);

      double relative_vinfall =  get_relative_infall_vel(block->xx[i_gal],block->xx[j_gal],
                                                        block->yy[i_gal],block->yy[j_gal],
                                                        block->zz[i_gal],block->zz[j_gal],
                                                        block->vxx[i_gal],block->vxx[j_gal],
                                                        block->vyy[i_gal],block->vyy[j_gal],
                                                        block->vzz[i_gal],block->vzz[j_gal],
                                                        distance);

      // Compute the quantities we need
      #pragma omp atomic
      block->n_pairs[i_bin] += 1;
      #pragma omp atomic
      block->sum_vinfall[i_bin] += relative_vinfall;
      #pragma omp atomic
      block->sum_vinfall_sq[i_bin] += (relative_vinfall*relative_vinfall);


      // If we considered all the neighbours of i_gal, it makes no sense to keep going
      // and we can jump to the next i_gal.
      ++n_igal_neighbours_so_far;
      if (n_igal_neighbours_so_far==block->n_neighbours[i_gal]) {
        if (params->verbose > 5)
          printf("We went through all the neighbours of i_gal=%lld, going to next i_gal.\n", i_gal);
        break;
      }  

    } // End of loop on j_gal

  } // End of loop on i_gal

  // Print to screen the total number of pairs
  if (params->verbose > 1)
    printf("Found %lld pairs that match the isolation criterion (%.4g percent of possible pairs).\n",
      n_pairs, 100.*n_pairs/((block->n_gals-1.)*block->n_gals/2.));

  // Compute the statistics for the block
  for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin) {

   // Average
    block->avg_vinfall[i_bin] = block->sum_vinfall[i_bin]/((double)(block->n_pairs[i_bin]));

    // The variance is <V_j^2> - <V_j>^2 , where j is a pair's index
    block->variance_vinfall[i_bin] = block->sum_vinfall_sq[i_bin]/((double)(block->n_pairs[i_bin]))
                                     - block->avg_vinfall[i_bin]*block->avg_vinfall[i_bin];

    if((block->variance_vinfall[i_bin] == 0.) && (block->n_pairs[i_bin]>1) && params->verbose>0)
      WARNINGPRINT("The bin %d for the block %d has zero variance but %lld pairs!\n", i_bin, block->id, block->n_pairs[i_bin]);

  }

  // Issue a warning if some point-like pairs were found
  if (params->verbose>0 && n_same_position_gals!=0)
    WARNINGPRINT("%lld point-like pairs of galaxies %s found and excluded!\n",
      n_same_position_gals, n_same_position_gals>1 ? "were":"was");

  return _SUCCESS_;

}
 

/* Read in the file associated to a block.
 * Here we allocate and fill the arrays block->xx, block->yy, block->zz,
 * block->vxx, block->vyy, block->vzz.
 * We also assign a value to block->n_gals.
 */
int data_block_read_file (
          params_struct *params,
          data_block_struct *block
          ) {


  // If 'params->n_gals_per_block=0', then we read all the galaxies in the block file.
  // Otherwise, we just read 'params->n_gals_per_block' galaxies per block.

  // Since 'read_six_column_csv_file' may overwrite the 'n_gals_read_from_file' argument,
  // we need to pass 'block->n_gals', not 'params->n_gals_per_block'
  block->n_gals_read_from_file = params->n_gals_per_block;

  // ==============================
  // = Read the file & apply cuts =
  // ==============================
  if (params->n_cuts==0) {
    zelda_call (read_six_column_csv_file (block->in_filename, &(block->n_gals_read_from_file),
                                  &(block->xx), &(block->yy), &(block->zz),
                                  &(block->vxx), &(block->vyy), &(block->vzz)));
    // Since there are no cuts, the number of galaxies considered
    // is equal to the number of galaxies read from the file.
    block->n_gals = block->n_gals_read_from_file;
  }
  else {
    // When applying the cuts, the variable 'block->n_gals' will contain just
    // the number of galaxies that survived the cuts
    zelda_call (read_and_cut_six_column_csv_file (block->in_filename, &(block->n_gals_read_from_file),
                                  &(block->xx), &(block->yy), &(block->zz),
                                  &(block->vxx), &(block->vyy), &(block->vzz),
                                  params->n_cuts,
                                  params->columns_to_cut,
                                  params->cuts,
                                  &(block->n_gals)
                                  ));
  }

  // Print some information
  if(block->n_gals_read_from_file==0) {
    printf("The block file is EMPTY: %s\n", block->in_filename);
  }
  else {
    if (params->verbose > 1) {    
      printf("Read %lld galaxies from file '%s'.\n", block->n_gals_read_from_file, block->in_filename);
      if (params->n_cuts > 0)
      printf("Of these, %lld survived the cuts (%.3g percent).\n", block->n_gals,
        100.*(double)(block->n_gals)/(double)(block->n_gals_read_from_file));    
    }
  }
  
  // Flag the block if zero galaxies survived the cuts
  if (block->n_gals==0) {
    block->is_empty_block = _TRUE_;
    return _SUCCESS_;
  }
   
   
  // ==========================
  // = Sorting (if requested) =
  // ========================== 
  if (params->sorting==_TRUE_){
    
    if(params->verbose > 0)
      printf("Sorting the block files according to their first column... ");
    
    // Get sorting indices for x column
    long long *sorted_indices;
    zelda_alloc(sorted_indices,sizeof(long long)*block->n_gals);
    if (get_sorting_indices(block->xx, block->n_gals, sorted_indices)==_FAILURE_) {
      ERRORPRINT("Something went wrong sorting the x coordinate\n","");
      return _FAILURE_;
    }
    
    // Do the actual sorting
    sort_array(block->xx,block->n_gals,sorted_indices);
    sort_array(block->yy,block->n_gals,sorted_indices);
    sort_array(block->zz,block->n_gals,sorted_indices);
    sort_array(block->vxx,block->n_gals,sorted_indices);
    sort_array(block->vyy,block->n_gals,sorted_indices);
    sort_array(block->vzz,block->n_gals,sorted_indices);
    
    if(params->verbose > 0)
      printf("done.\n");
  }

  return _SUCCESS_;

}



/* Save the statistic associated to a block to its output file */
int data_block_write_file (
          params_struct *params,
          data_block_struct *block
          ) {

  int i_bin=0;

  FILE *out_file = fopen(block->out_filename,"w");
  if (out_file==NULL) {
    printf("ERROR, %s: Cannot open file %s for writing.\n", __func__, block->out_filename);
    return _FAILURE_;
  }

  if (params->verbose > 1)
    printf("Storing results for block %d in file %s... ", block->id+1, block->out_filename);

  char labels[][128] = {
    "edges",
    "midpts",
    "npairs",
    "sumV",
    "avgV",
    "sumV^2",
    "varV",
    "sigmaV"
  };

  fprintf(out_file, "# %28s %30s %30s %30s %30s %30s %30s %30s\n",
    labels[0], labels[1], labels[2], labels[3], labels[4], labels[5], labels[6], labels[7]);

  for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin) {

    fprintf(out_file, "%+30f %+30f %30lld %+30f %+30f %+30f %+30f %+30f\n", 
      params->r_edges[i_bin],
      params->r_midpoints[i_bin],
      block->n_pairs[i_bin],
      block->sum_vinfall[i_bin],
      block->avg_vinfall[i_bin],
      block->sum_vinfall_sq[i_bin],
      block->variance_vinfall[i_bin],
      sqrt(block->variance_vinfall[i_bin])
      );

  }

  if (params->verbose > 1)
    printf("Done.\n"); fflush(stdout);
  
  fclose(out_file);
  
  return _SUCCESS_;

}



/* Determine the number of neighbours for each galaxy in a sphere of radius
 * 'r_isolation' centred on the galaxy.
 * In the worst case (big 'r_isolation' and unlimited 'n_neighbours_max'),
 * this is a n_gals^2 algorithm.  However, for reasonable values of the
 * parameters, it is much faster.  This is possible because the galaxies
 * are sorted according to their 'x' component.
 */
int data_block_neighbours (
        params_struct *params,
        data_block_struct *block
        )
{

  long long i_gal=0, j_gal=0, gals_on_edges=0;

  // 'n_neighbours' will contain for each galaxy the number of galaxies in a sphere
  // or radius params->r_isolation centered in that galaxy
  zelda_calloc(block->n_neighbours, block->n_gals, sizeof(long long));

  // Handy shortcut
  double r_iso = params->r_isolation;

  if (params->verbose > 3)
    printf("%s: Entering the main loop.\n", __func__);
    
  // ==========================
  // = Main cycle on galaxies =
  // ==========================
  //  Beginning of parallel region
  # pragma omp parallel for private(i_gal, j_gal)
  
  // We shall exploit the fact that the galaxies have been sorted in the x-direction
  // to allow a faster computation.  In order to do so, we split the loop on the j_gals
  // in a forward and a backward loop.

  // *** Loop in the forward direction ***
  for (i_gal=0; i_gal < block->n_gals; ++i_gal) {        

    // Condition to remove edge effects. That is, we don't want to have pairs that are
    // in the edges of boxes because it may be that they are not really
    // isolated pairs but have a nearby galaxy on the block next to it. 
    // This condition has nothing to do with the number of neighbours, but it is convenient
    // to put it here, since removing galaxies at this early stage saves some computational time.
    // N.B.  The numbering of i_x, i_y, i_z starts from 1, hence the -1's below.
    if (params->remove_galaxies_on_the_edges == _TRUE_) {
      if(  (block->xx[i_gal] < (block->i_x-1)*block->side_length + r_iso) || (block->xx[i_gal] > block->i_x*block->side_length - r_iso) 
        || (block->yy[i_gal] < (block->i_y-1)*block->side_length + r_iso) || (block->yy[i_gal] > block->i_y*block->side_length - r_iso)
        || (block->zz[i_gal] < (block->i_z-1)*block->side_length + r_iso) || (block->zz[i_gal] > block->i_z*block->side_length - r_iso))
    {
        
        block->n_neighbours[i_gal] = _BADFLAG_;
    
        #pragma omp atomic
        ++gals_on_edges;
    
        continue;
      }
    }

  
  // *** Loop in the backward direction ***
    // IMPORTANT:  YOU MAY CONSIDER REMOVING THE BACKWARD LOOP,  AND USING SOMETHING LIKE THIS IN
    // THE BODY OF THE FORWARD LOOP:
    // if (distance < r_iso) {
    //   ++block->n_neighbours[i_gal];
    //   ++block->n_neighbours[j_gal];      
    // }  

    // Note that if i_gal = 0 (first iteration) we are not entering this loop    
    for (j_gal=i_gal-1; j_gal>=0; --j_gal) {
    
      // Distance in x between the pair
      double delta_x = block->xx[i_gal]-block->xx[j_gal];
  
      // Check that x is sorted
      if (delta_x < 0)
        WARNINGPRINT("The input file %s has not been sorted as it should.\n",block->in_filename);
    
      // Stop looking for pairs when the delta_x is bigger than the isolation criteria.
      // (Remember that 'block->n_neighbours[i_gal]' is zero by default.)
      if( delta_x > r_iso )
        break;
  
      double distance = get_distance(block->xx[i_gal],block->xx[j_gal],
                                  block->yy[i_gal],block->yy[j_gal],
                                  block->zz[i_gal],block->zz[j_gal]);
  
      // If j_gal is within a 'r_isolation' sphere from 'i_gal', then 
      // the two galaxies are neighbours.
      if (distance < r_iso)
        ++block->n_neighbours[i_gal];
  
      // Break if the sphere becomes too crowded, and flag the galaxy as a bad one.
      if (block->n_neighbours[i_gal] > params->n_neighbours_max ) {
        block->n_neighbours[i_gal] = _BADFLAG_;
        break;
      }
  
    } // End on backward loop
  
    // If the backward semi-sphere is already too crowded, then it makes
    // no sense to go on
    if (block->n_neighbours[i_gal] == _BADFLAG_)
      continue;
  
  
    // =================================
    // = Loop in the forward direction =
    // =================================
  
    // Note that if i_gal = block->n_gals-1 (last iteration) we are not entering this loop
    for (j_gal=i_gal+1; j_gal<block->n_gals; ++j_gal) {
  
      // Distance in x between the pair
      double delta_x = block->xx[j_gal]-block->xx[i_gal];
  
      // Check that x is sorted
      if (delta_x < 0)
        WARNINGPRINT("The input file %s has not been sorted as it should.\n",block->in_filename);
  
      // Stop looking for pairs when the delta_x is bigger than the isolation criterion.
      // (Remember that 'block->n_neighbours[i_gal]' is zero by default.)      
      if( delta_x > r_iso )
        break;
  
      double distance = get_distance(block->xx[i_gal],block->xx[j_gal],
                                  block->yy[i_gal],block->yy[j_gal],
                                  block->zz[i_gal],block->zz[j_gal]);
  
      // n_neighbours[i_gal] is in the range (0,n_neighbours_max) from the previous backward loop.
      // If it were bigger than 'n_neighbours_max', then we wouldn't be inside this loop.
      if (distance < r_iso)
        ++block->n_neighbours[i_gal];
  
      // Break if the sphere becomes too crowded
      if (block->n_neighbours[i_gal] > params->n_neighbours_max) {
        block->n_neighbours[i_gal] = _BADFLAG_;
        break;
      }
  
    } // End of the forward loop on j_gal
    
  } // End of loop on i_gal

  if ((params->verbose > 1) && (params->remove_galaxies_on_the_edges==_TRUE_)) {
    printf("Expected percentage of galaxies near the edges = %g\n", (1-pow(1 - 2*r_iso/block->side_length,3))*100);
    printf("The percentage of galaxies near the edges is = %g \n", (double)gals_on_edges/(double)(block->n_gals)*100.);
  }

  // Debug.  Print the number of neighbours of those galaxies wih less
  // than 'n_neighbours_max' neighbours.
  // for(i_gal=0; i_gal<block->n_gals; ++i_gal) {
  //   if(block->n_neighbours[i_gal]!=_BADFLAG_)
  //     printf("block->n_neighbours[%lld-1]=%lld\n", i_gal+1, block->n_neighbours[i_gal]);
  // }

  return _SUCCESS_;
  
}


/* Free temporary arrays, i.e. all those that we do not want to keep
 * after we finished processing the block.
 */
int data_block_free_workspace (
          data_block_struct *block
          ) {

  free(block->xx);
  free(block->yy);
  free(block->zz);
  free(block->vxx);
  free(block->vyy);
  free(block->vzz);
  free(block->n_neighbours);

  return _SUCCESS_;

}



/* Free the memory associated to the structure */
int data_block_free (
        data_block_struct *block
        ) {
          
  free(block->n_pairs);
  free(block->sum_vinfall);
  free(block->avg_vinfall);
  free(block->sum_vinfall_sq);
  free(block->variance_vinfall);

  return _SUCCESS_;

}



/* Calculate pairwise velocity statistics considering only isolated pairs of galaxies.
 * The results obtained with this function should be equal to those obtained with
 * 'data_block_vinfall_stat_with_isolation' when params->n_neighbours_max=1.
 * The difference is that here we use a different algorithm, useful for testing purposes.
 * The speed of the two algorithms is more or less the same.
 */
int data_block_vinfall_stat_with_total_isolation (
          params_struct *params,
          data_block_struct *block
          ) {


  long long i_gal=0, j_gal=0;
  int i_bin=0;

  // Handy shortcut
  double r_iso = params->r_isolation;

  // This counter will contain for each galaxy the number of galaxies in a sphere
  // or radius r_iso centered in that galaxy
  long long n_pairs_in_sphere = 0;
  long long *neighbours_indxs;
  zelda_alloc(neighbours_indxs, sizeof(long long)*block->n_gals);

  if (params->verbose > 3)
    printf("%s: Entering the main loop.\n", __func__);
    
  long long gals_on_edges = 0;
  
  // ==========================
  // = Main cycle on galaxies =
  // ==========================
  //  Beginning of parallel region
  # pragma omp parallel for private(n_pairs_in_sphere, i_gal, j_gal)
  
  for (i_gal=0; i_gal < block->n_gals; ++i_gal) {

    // Condition to remove edge effects. That is, we don't want to have pairs that are
    // in the edges of boxes because it may be that they are not really
    // isolated pairs but have a nearby galaxy on the block next to it. 
    if (params->remove_galaxies_on_the_edges == _TRUE_) {
      if(  (block->xx[i_gal] < (block->i_x-1)*block->side_length + r_iso) || (block->xx[i_gal] > block->i_x*block->side_length - r_iso) 
        || (block->yy[i_gal] < (block->i_y-1)*block->side_length + r_iso) || (block->yy[i_gal] > block->i_y*block->side_length - r_iso)
        || (block->zz[i_gal] < (block->i_z-1)*block->side_length + r_iso) || (block->zz[i_gal] > block->i_z*block->side_length - r_iso))
      {

        #pragma omp atomic
        ++gals_on_edges;

        neighbours_indxs[i_gal] = _BADFLAG_;
        
        continue;
      }
    }
     
             
    // Reset the counter
    n_pairs_in_sphere = 0;

    // ==================================
    // = Loop in the backward direction =
    // ==================================
    // Note that if i_gal = 0 (first iteration) we are not entering this loop    
    for (j_gal=i_gal-1; j_gal>=0; --j_gal) {
    
      // Distance in x between the pair
      double delta_x = block->xx[i_gal]-block->xx[j_gal];

      // Check that x is sorted
      if (delta_x < 0)
        WARNINGPRINT("The input file %s has not been sorted as it should.\n",block->in_filename);
    
      // Stop looking for pairs when the delta_x is bigger than the isolation criteria
      if( delta_x > r_iso )
        break;

      double distance = get_distance(block->xx[i_gal],block->xx[j_gal],
                                  block->yy[i_gal],block->yy[j_gal],
                                  block->zz[i_gal],block->zz[j_gal]);

      // Enforce the isolation criterion
      if (distance < r_iso) {
        ++n_pairs_in_sphere;
        neighbours_indxs[i_gal] = j_gal;
      }

      // Break if the sphere becomes too crowded
      if (n_pairs_in_sphere > 1) {
        neighbours_indxs[i_gal] = _BADFLAG_;
        break;
      }

    } // End on backward loop
  
    // If the backward semi-sphere is already too crowded, then it makes
    // no sense to go on
    if (neighbours_indxs[i_gal] == _BADFLAG_)
      continue;


    // =================================
    // = Loop in the forward direction =
    // =================================

    // Note that if i_gal = block->n_gals-1 (last iteration) we are not entering this loop
    for (j_gal=i_gal+1; j_gal<block->n_gals; ++j_gal) {

      // Distance in x between the pair
      double delta_x = block->xx[j_gal]-block->xx[i_gal];

      // Check that x is sorted
      if (delta_x < 0)
        WARNINGPRINT("The input file %s has not been sorted as it should.\n",block->in_filename);

      // Stop looking for pairs when the delta_x is bigger than the isolation criterion
      if( delta_x > r_iso )
        break;
      
      double distance = get_distance(block->xx[i_gal],block->xx[j_gal],
                                  block->yy[i_gal],block->yy[j_gal],
                                  block->zz[i_gal],block->zz[j_gal]);

      // n_pairs_in_sphere is either 0 or 1 from the previous backward loop
      // (if it were 2, then we wouldn't be inside this loop.)
      if (distance < r_iso) {
        ++n_pairs_in_sphere;
        neighbours_indxs[i_gal] = j_gal;
      }

      // Break if the sphere becomes too crowded
      if (n_pairs_in_sphere > 1) {
        neighbours_indxs[i_gal] = _BADFLAG_;
        break;
      }

    } // End of the forward loop on j_gal

    // Account for the case where no neighbouring galaxies are found
    if (n_pairs_in_sphere==0)
      neighbours_indxs[i_gal] = _BADFLAG_;
  
  } // End of loop on i_gal
  


  // ======================
  // = Get the statistics =
  // ======================
  // Allocate memory for the result arrays
  block->sum_vinfall = (double*)calloc(params->n_separation_bins,sizeof(double));
  block->avg_vinfall = (double*)calloc(params->n_separation_bins,sizeof(double));
  block->sum_vinfall_sq = (double*)calloc(params->n_separation_bins,sizeof(double));
  block->variance_vinfall = (double*)calloc(params->n_separation_bins,sizeof(double));

  // This variable may be useful if we are dealing with log spaced separation bins.
  // We compute it outside the main loop.
  double log_r_min;
  if (params->binning_mode == LOG_BINNING)
    log_r_min = log(params->r_min);
  
  // Accumulator to count the pairs with zero distance
  long long n_same_position_gals=0;
  
  // Process the neighbours_indxs array to find actual pair that satisfy the criterion,
  // and store the statistics
  for(i_gal=0; i_gal<block->n_gals; ++i_gal) {

    // =================
    // = Find the pair =
    // =================

    // If i_gal has zero or too many neighbours, we are not interested in collecting
    // any statistics
    if (neighbours_indxs[i_gal] == _BADFLAG_)
      continue;

    // We already know that the only neighbour of i_gal is j_gal ...
    j_gal = neighbours_indxs[i_gal];
    
    // ... so we only need to check if the only neighbour of j_gal is i_gal!
    if (neighbours_indxs[j_gal] != i_gal)
      continue;

    // When we get here, we are sure that i_gal and j_gal form a pair. However, for each pair
    // we'll get here twice because the loop is on i_gal=0:n_gals.  Therefore, we jump to the
    // next iteration when j_gal is smaller than i_gal.
    if(j_gal < i_gal)
      continue;
      
    // =====================
    // = Do the statistics =
    // =====================
    double distance = get_distance(block->xx[i_gal],block->xx[j_gal],
                                  block->yy[i_gal],block->yy[j_gal],
                                  block->zz[i_gal],block->zz[j_gal]);

    // Perform a check on point-like pairs, i.e. repeated lines in the input files
    if(distance==0.) {
      if (params->verbose>3) {
        WARNINGPRINT("Pair (%lld,%lld) is point-like!\n",i_gal+1,j_gal+1);
        if( (block->vxx[i_gal]!=block->vxx[j_gal]) )
          printf("         (the velocity is different though!)\n");
      }
      #pragma omp atomic
      ++n_same_position_gals;
      continue;
    }


    // Check that the pair is in the requested range, otherwise skip it.
    if ((distance < params->r_edges[0]) || (distance >= params->r_edges[params->n_separation_bins]))
      continue;
      
    // Assign the pair to its separation bin.
    if (params->binning_mode == LIN_BINNING)
      i_bin = (int)((distance-params->r_min)/params->r_bin_width);
      
    else if (params->binning_mode == LOG_BINNING)
      i_bin = (int)((log(distance) - log_r_min)/params->r_bin_width);
      
    else if (params->binning_mode == CUSTOM_BINNING) {
      for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin)
        if ( (distance >= params->r_edges[i_bin]) && (distance < params->r_edges[i_bin+1]) )
          continue;
    }

                                  
    double relative_vinfall =  get_relative_infall_vel(block->xx[i_gal],block->xx[j_gal],
                                                      block->yy[i_gal],block->yy[j_gal],
                                                      block->zz[i_gal],block->zz[j_gal],
                                                      block->vxx[i_gal],block->vxx[j_gal],
                                                      block->vyy[i_gal],block->vyy[j_gal],
                                                      block->vzz[i_gal],block->vzz[j_gal],
                                                      distance);

    // Compute the quantities we need
    #pragma omp atomic    
    block->n_pairs[i_bin] += 1;
    #pragma omp atomic
    block->sum_vinfall[i_bin] += relative_vinfall;
    #pragma omp atomic
    block->sum_vinfall_sq[i_bin] += (relative_vinfall*relative_vinfall);

    if (params->verbose > 3) {
      printf("Galaxies (%lld,%lld)-1 are matched.\n", i_gal+1, j_gal+1);
      printf("    x=%20g, y=%20g, z=%20g\n", block->xx[i_gal], block->yy[i_gal], block->zz[i_gal]);
      printf("    x=%20g, y=%20g, z=%20g\n", block->xx[j_gal], block->yy[j_gal], block->zz[j_gal]);        
    }


  } // End of loop on the elements of neighbours_indxs

  // Print to screen the total number of pairs
  if (params->verbose > 1) {
    long long n_total_pairs=0;
    for(i_bin=0; i_bin<params->n_separation_bins; ++i_bin)
      n_total_pairs += block->n_pairs[i_bin];
    printf("Found %lld pairs that match the isolation criterion (%.4g percent of possible pairs).\n",
      n_total_pairs, 100.*n_total_pairs/((block->n_gals-1.)*block->n_gals/2.)) ;
  }
  
  if ((params->verbose > 1) && (params->remove_galaxies_on_the_edges==_TRUE_)) {
    printf("Expected percentage of galaxies near the edges = %g\n", (1-pow(1 - 2*r_iso/block->side_length,3))*100);
    printf("The percentage of galaxies near the edges is = %g \n", (double)gals_on_edges/(double)(block->n_gals)*100.);
  }


  // Compute the statistics for the block
  for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin) {
  
   // Average
    block->avg_vinfall[i_bin] = block->sum_vinfall[i_bin]/((double)(block->n_pairs[i_bin]));
  
    // The variance is <V_j^2> - <V_j>^2 , where j is a pair's index
    block->variance_vinfall[i_bin] = block->sum_vinfall_sq[i_bin]/((double)(block->n_pairs[i_bin]))
                                     - block->avg_vinfall[i_bin]*block->avg_vinfall[i_bin];
   
    if((block->variance_vinfall[i_bin] == 0.) && (block->n_pairs[i_bin]>1) && params->verbose>0)
      WARNINGPRINT("The bin %d for the block %d has zero variance but %lld pairs!\n", i_bin, block->id, block->n_pairs[i_bin]);
    
  }
  
  // We are done!
  free(neighbours_indxs);
  
  // Issue a warning if some point-like pairs were found
  if (params->verbose>0 && n_same_position_gals!=0)
    WARNINGPRINT("%lld point-like pairs of galaxies %s found and excluded!\n",
      n_same_position_gals, n_same_position_gals>1 ? "were":"was");

  return _SUCCESS_;

}