#ifndef _DATABLOCKTOOLS_
#define _DATABLOCKTOOLS_

#include "params.h"

/* Structure containing all that is needed to compute the statistics
 * of a single block of data.
 */
typedef struct {
  
  // Block identifier  
  int id;
  
  // Length of the side of the block. Used to remove the edges effects.
  double side_length;
  
  // Bins in x,y,z corresponding to the block
  int i_x;
  int i_y;
  int i_z;    
  
  // Input file associated to the block
  FileName in_filename;

  // Output file associated to the block
  FileName out_filename;  
    
  // Number of galaxies considered for the data block.  This number
  // is equal to the number of galaxies read from the input file,
  // minus the galaxies that did not survive the cuts
  long long n_gals;

  // Number of galaxies read from the input file.  Consider that if 'params->n_gals_per_block'
  // is set to something different from zero, then only that amount of
  // galaxies will be read from the file (useful for quick debug).
  long long n_gals_read_from_file;
    
  // Array with the number of pairs in each bin (result array)
  long long *n_pairs;

  // Array that accumulates the relative infall velocities of pairs of galaxies
  // for each bin (result array)  
  double *sum_vinfall;

  // Average of the relative infall velocities for each bin. It is just
  // sum_vinfall divided by the number of pairs in each bin (result array)
  double *avg_vinfall;
  
  // Array needed to compute the total variance of the dataset
  double *sum_vinfall_sq;
  
  // Variance of v_infall for the single block
  double *variance_vinfall;
  
  // Number of neighbours for each galaxy.  This array is used to enforce the
  // isolation criterion with an arbitrary number of neighbours, and is filled
  // by the function 'data_block_neighbours'
  long long *n_neighbours;

  // Flag that identifies the block has empty.  This happens if no galaxy in
  // in the block survived the cuts made in 'data_block_read_file'.
  int is_empty_block;
  
  // Position and velocities of galaxies in the block
  // These are temporary/workspace variables, in the sense that they are freed
  // after the block has been processed
  double *xx;
  double *yy;
  double *zz;
  double *vxx;
  double *vyy;
  double *vzz;
  
      
    
} data_block_struct;



/*  Fill a data_block_struct structure with the data contained in block file,
 *  and do all is needed to do on it
 */
int data_block_init (
        params_struct *params,
        data_block_struct *block,
        int id
        );
                

/*  Read the file corresponding to the block */
int data_block_read_file (
        params_struct *params,
        data_block_struct *block
        );

/* Write the results obtained by processing a block into its file */
int data_block_write_file (
        params_struct *params,
        data_block_struct *block
        );

        

/*  Compute the pairwise velocity statistics for all the pairs in the block */
int data_block_vinfall_stat (
        params_struct *params,
        data_block_struct *block
        );


/* Compute the pairwise velocity statistics for isolated galaxy pairs.
 * The degree of isolation is controlled by 'params->r_isolation' and
 * 'params->n_neighbours_max'  */
int data_block_vinfall_stat_with_isolation (
        params_struct *params,
        data_block_struct *block
        );

/* Same as above, but selecting only those pairs that are completely isolated (i.e. 0 neighbours) */
int data_block_vinfall_stat_with_total_isolation (
        params_struct *params,
        data_block_struct *block
        );


/* Determine the number of neighbours for each galaxy */
int data_block_neighbours (
        params_struct *params,
        data_block_struct *block
        );

/* Free temporary arrays */
int data_block_free_workspace (
        data_block_struct *block
        );


/* Free the memory associated to the structure */
int data_block_free (
        data_block_struct *block
        );





#endif