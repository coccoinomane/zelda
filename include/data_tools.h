#ifndef _DATATOOLS_
#define _DATATOOLS_

#include "data_block_tools.h"

/* Data structure.  It contains all the blocks and the statistic associated
 * to the whole dataset
 */
typedef struct {
  
  // Blocks of galaxies
  int n_blocks;
  data_block_struct *blocks;

  // Filename where to write the results
  char out_filename[NAME_MAX];

  // Total number of galaxies
  long long n_gals;


  // =========================
  // = Statistical variables =
  // =========================
  // Array with the number of pairs in each bin
  long long *n_pairs;
  long long n_pairs_total;


  /* Pair-averaged quantities */
  // These quantities are related to the total number of pairs in the dataset,
  // in the sense that the sums & expection values are taken over that ensemble.
  // This means that the random variable is the infall velocity of each pair.
  // More details in the implementation of the function data_vinfall_stats.

  // Arrays that accumulates the relative infall velocities (and their squared
  // versions) of pairs of galaxies for each bin
  double *sum_vinfall;
  double *sum_vinfall_sq;
  
  // Average of the relative infall velocities for each bin. It is just
  // sum_vinfall divided by the number of pairs in each bin
  double *avg_vinfall;

  // Variance of v_infall for the whole dataset
  double *variance_vinfall;

  /* Block-averaged quantities */
  // These quantities are related to the number of blocks in the dataset,
  // in the sense that the sums & expection values are taken over that ensemble.
  // This means that the random variable is the (mean) infall velocity of each block.
  // More details in the implementation of the function data_vinfall_stats.
  // The meaning of each of these quantities is equivalent to those of the 
  // pair averaged quantities.
  double *sum_vinfall_means;
  double *sum_vinfall_means_sq;
  double *avg_vinfall_means;
  double *variance_vinfall_means;  
  double *averaged_variance_vinfall;
  

} data_struct;



/* Initialise a data structure
 */
int data_init (
        params_struct *params,
        data_struct *data
        );


/* Compute the infall velocities statistics related to full dataset,
 * starting from those already calculated for the single blocks
 */
int data_vinfall_stat (
        params_struct *params,
        data_struct *data
        );

/* Write to file the results for the whole dataset */
int data_write_file (
          params_struct *params,
          data_struct *data
          );



/* Read in the pre-computed files containing the linear infall velocities */
int data_read_vinfall_stat_from_files (
        params_struct *params,
        data_struct *data
        );

/* Free the memory associated to the structure */
int data_free (
        data_struct *data
        );


#endif