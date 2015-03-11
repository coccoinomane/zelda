/** test_cuts.c
 * Created by Guido W. Pettinari on 20/11/2011
 * Last modified by Guido W. Pettinari on 20/11/2011
 *
 * Test the function 'read_and_cut_six_column_csv_file' on the file
 * 'tests/test_cut_id_type_np_mvir_haloid_rmag.dat.txt'.
 * With three arguments, apply only one cut.  The first argument
 * is the column with respect to which the cut should be applied.
 * The second and third arguments are the minimum and maximum allowed
 * value for that column, respectively.
 * With seven arguments, apply two cuts. 
 * As an example, try these commands:
 *
 * ./test_cuts 3 40 60
 *
 * Output:
 * n_lines_in_file = 20
 * n_lines_kept = 3
 * 5557.000000,2.000000,53.000000,4.131152,0.000000,-17.128418
 * 5824.000000,2.000000,49.000000,4.991809,0.000000,-17.368528
 * 6076.000000,2.000000,60.000000,6.196728,0.000000,-17.736202
 *
 * ./test_cuts 3 40 60 1 5000 6000
 *
 * Output:
 * n_lines_in_file = 20
 * n_lines_kept = 2
 * 5557.000000,2.000000,53.000000,4.131152,0.000000,-17.128418
 * 5824.000000,2.000000,49.000000,4.991809,0.000000,-17.368528
 *
 */
 
#include "tools.h"
#include "stdio.h"
#include "stdlib.h"

int main (int argc, char const *argv[]) {

  // Parse arguments
  if ((argc != 4) && (argc != 7)) {
    printf("usage:     %s <column_idx> <min> <max> [<column_idx> <min> <max>]\n", argv[0]);
    return 1;
  }

  // File to be read
  char file_name[256] = "tests/test_cut_id_type_np_mvir_haloid_rmag.dat.txt";

  // Cut-related arguments for the function 'read_and_cut_six_column_csv_file'
  int n_cuts;
  int * which_columns_to_cut;
  double * cuts;
  
  // Cycle variable
  int i_cut=0;
  
  if (argc == 4)
    n_cuts = 1;
  else
    n_cuts = 2;
    
  cuts = malloc(2*n_cuts*sizeof(double));
  which_columns_to_cut = malloc(sizeof(long long)*n_cuts);
  which_columns_to_cut[0] = atoi(argv[1]);
  cuts[0] = atof(argv[2]);
  cuts[1] = atof(argv[3]);
  
  // If 2 cuts are requested, fill the rest of the arrays
  if (n_cuts == 2) {
    which_columns_to_cut[1] = atoi(argv[4]);
    cuts[2] = atof(argv[5]);
    cuts[3] = atof(argv[6]);    
  }

  // Print some info
  for(i_cut=0; i_cut<n_cuts; ++i_cut) {
    printf("\nWe shall perform the following cuts:\n");
    printf("    column %d between %g and %g\n", which_columns_to_cut[i_cut], cuts[2*i_cut], cuts[2*i_cut+1]);
  } 
  printf("\n");
      
  // Output variables
  long long n_lines_in_file;
  long long n_lines_kept;  
  double *c1, *c2, *c3, *c4, *c5, *c6;

  read_and_cut_six_column_csv_file (
          file_name,                                         // In
          &n_lines_in_file,                                  // In & out
          &c1, &c2, &c3, &c4, &c5, &c6,                      // Out
          n_cuts,                                            // In
          which_columns_to_cut,                              // In
          cuts,                                              // In
          &n_lines_kept                                      // Out
          );  

  printf("n_lines_in_file = %lld\n", n_lines_in_file);
  printf("n_lines_kept = %lld\n", n_lines_kept);
  long long i=0;
  for (i=0; i<n_lines_kept; i++) {
    printf("%lf,%lf,%lf,%lf,%lf,%lf\n", c1[i], c2[i], c3[i], c4[i], c5[i], c6[i]);
  }
  
  return 0;
  
}