#include "zelda.h"

int main(int argc, char ** argv){
  

  if (argc != 5) {
     printf("usage:  %s <input_filename> <format for sprintf> <size_of_grid> <n_bins_grid>\n", argv[0]);
     return _FAILURE_;
  }
  char * in_filename = argv[1];
  char * format = argv[2];
  double size_grid = atof(argv[3]);
  int n_bins_grid = atoi(argv[4]);
  
  if(split_six_column_csv_file (in_filename, format, size_grid, n_bins_grid) == _FAILURE_) {
    printf("Something went wrong reading and splitting the file\n");
    return _FAILURE_;
  }
  
  int n_out_files = n_bins_grid*n_bins_grid*n_bins_grid; 
  printf("Congratulations! The file %s is now split into %d files\n", in_filename, n_out_files);
   
  return _SUCCESS_;

}
