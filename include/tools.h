/** tools.h
 * Last modified by Guido W. Pettinari on 20/11/2011
 *
 * Support functions for the Zelda code.
 *
 */

#ifndef _TOOLS_
#define _TOOLS_


/* Build a format specifier for sprintf-like functions, so that only the
 * n-th element of a string is considered.  Elements are assumed to be 
 * delimited by 'delimiter'.
 * An example is the easiest way to illustrate this function.
 * entry_type = 'd'
 * delimiter = ','
 * n = 3
 * -> out = "%*d,%*d,%d
 */
int format_for_nth_column (
       char entry_type[],
       char delimiter[],
       int n,
       char * out
       );


/* Given a comma-separated-value text file containing six or more entries per line,
 * apply the specified cuts and store the surviving lines of the first six columns
 * in the six specified arrays.
 *
 * The amount of lines read in is given by the input/output argument n_lines_pointer.
 * If it is zero, read all the lines in the file until EOF is found.
 *
 * The cuts should be specified using the following variables (see below for details)
 * int n_cuts
 * int *which_columns_to_cut
 * double ** cuts
 */
int read_and_cut_six_column_csv_file (
        char *file_name,                                             // In
        long long *n_lines_read_from_file,                           // In & out
        // pc1, pc2, ... are pointers to the 6 arrays that will
        // contain the data. Notation: pcn = pointer to column 'n'  
        double **pc1, double **pc2, double **pc3,                    // Out
        double **pc4, double **pc5, double **pc6,                    // Out
        // Number of cuts to apply
        int n_cuts,                                                  // In
        // Index of the columns to be used for the cuts
        int * which_columns_to_cut,                                  // In
        // Range of the cuts:
        // min[i_cut] = cuts[2*i_cut]
        // max[i_cut] = cuts[2*i_cut+1]
        double * cuts,                                              // In
        // Number of lines kept after the cuts. It is also the
        // size of the pcn vectors.
        long long *n_lines_kept                                      // Out
        );






/* Given a comma-separated-value text file containing six or more entries per line,
 * store the first six columns in the six specified arrays.
 * The amount of lines read in is given by the input/output argument n_lines_pointer.
 * If it is zero, read all the lines in the file until EOF is found.
 */
int read_six_column_csv_file (
        // IN
        char *file_name,
        // OUT
        long long *n_lines_pointer,
        // Will contain the data in arrays
        // Notation: pcn = pointer to column 'n'
        double **pc1, double **pc2, double **pc3,
        double **pc4, double **pc5, double **pc6
        );






/* Given a 3+ column file, split/bin it into blocks according to the values of
 * the first three columns.  We us a linear grid with origin in the point (0,0,0).
 */
int split_six_column_csv_file (
        char *in_filename,
        char *format,
        double size_grid,
        int n_bins_grid
        );






// Indexes an array array[0..n-1], i.e., outputs the array indx[0..n-1] such that array[indx[j]] 
// is in ascending order for j = 0, 1, 2, . . . , n-1 by the use of a quicksort algorithm. 
// The input quantities n and array are not changed.
int get_sorting_indices(double *array, long long n, long long *indx);

/* Function that takes in an array and an array of sorting indices and outputs the sorted array
*/
int sort_array(double *array, long long n, long long *sorting_indices);


// Given two positions, get the distance between them
double get_distance(
  		  //IN
        double x1, double x2, double y1, double y2, double z1, double z2 );

// Given the positions and velocities of two objects, get the relative infall velocity between them
double get_relative_infall_vel(
        //IN
        double x1, double x2, double y1, double y2, double z1, double z2,
        double vx1, double vx2, double vy1, double vy2, double vz1, double vz2,
        double distance);
      
      
// Return an array with linearly spaced points
double * lin_space (double x_min, double x_max, int n_points);

// Return an array of logarithmically spaced points        
double * log_space (double x_min, double x_max, int n_points);    




#endif