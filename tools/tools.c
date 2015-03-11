#include "tools.h"
#include "common.h"


// =========================
// = format_for_nth_column =
// =========================
int format_for_nth_column (
        char entry_type[],
        char delimiter[],
        int n,
        char * out
        ) {
         
  int i=0;

  // Check the input
  zelda_test (n<=0, "the column number 'n' should be greater than zero.")
  
  // The first n-1 columns need to be ignored.  To do so, we use the
  // assignment-suppression character *.
  char ignore[128];
  strcpy(ignore, "%*");
  strcat(ignore, entry_type);
  strcat(ignore, delimiter);  
  
  // We shall get only the n-th column
  char get[128];
  strcpy(get, "%");
  strcat(get, entry_type);

  // If we want the first column, then we do not need to ignore anything
  if (n==1) {
    strcpy(out, get);
    return _SUCCESS_;
  }

  // Build the format specifier
  strcpy(out, ignore);
  for(i=1; i<n-1; ++i) {
    strcat(out, ignore);
  }
  strcat(out, get);
  
  // Debug
  // printf("get = %s\n", get);
  // printf("ignore = %s\n", ignore);
  // printf("out = %s\n", out);  
  
  return _SUCCESS_;
}


/// ====================================
/// = read_and_cut_six_column_csv_file =
/// ====================================
int read_and_cut_six_column_csv_file (
        char *file_name,                                                  // In
        long long *n_lines_read_from_file,                                // In & out
        // pc1, pc2, ... are pointers to the 6 arrays that will
        // contain the data. Notation: pcn = pointer to column 'n'  
        double **pc1, double **pc2, double **pc3,                         // Out
        double **pc4, double **pc5, double **pc6,                         // Out
        // Cuts to apply
        int n_cuts,                                                       // In
        int * which_columns_to_cut,                                       // In
        double * cuts,                                                    // In
        long long *n_lines_kept                                           // Out
        ) {
    
  // Cycle variables
	long long i_line=0;
  int i_cut=0;

	// Open the data file
	FILE *fp = fopen(file_name,"r");
	if (fp == NULL) {
    ERRORPRINT("Couldn't open file %s.\n", file_name);
    return _FAILURE_;
  }
	  
  // Buffer that will contain a line of the file
  char line[LINE_MAX];

	// =========================================
	// = Count the number of lines in the file =
	// =========================================
  // If the user gave a number different from 0 as n_lines, then she wants only to 
  // consider the lines up to that number.  Otherwise, she wants all the lines in
  // the file, and we need to count them.
  if (*n_lines_read_from_file==0) {
  	while ( fgets(line, LINE_MAX, fp) != NULL )
  			++(*n_lines_read_from_file);
  }


  // ==================
  // = Apply the cuts =
  // ==================
  // Rewind the file
  rewind(fp);

  // When we will apply the cuts, we shall remember which lines should be kept
  // using the logical array 'keep_line'
  short * keep_line;
  zelda_alloc (keep_line, (*n_lines_read_from_file)*sizeof(short));

  // Create the format specifiers to select the columns we are interested into.
  // We assume that all the entries in the file are in the form of a double
  // precision number.
  char format[n_cuts][512];

  // printf("\nWe shall perform the following cuts:\n");
  for(i_cut=0; i_cut<n_cuts; ++i_cut) {
    // printf("    column %d between %g and %g\n", which_columns_to_cut[i_cut], cuts[2*i_cut], cuts[2*i_cut+1]);
    zelda_call(format_for_nth_column ("lf", ",", which_columns_to_cut[i_cut], format[i_cut]));
  }
  // printf("\n");

  // Initialise the number of kept lines to zero
  *n_lines_kept = 0;
  
  // **** Loop on the lines & on the cuts ****
  for (i_line=0; i_line<*n_lines_read_from_file; i_line++) {
  	
  	// Read a line
    char *fgets_output = fgets(line, LINE_MAX, fp);

    // By default, assume that the line will survive the cuts
    keep_line[i_line] = _TRUE_;

    // Apply all the cuts sequentially
    for(i_cut=0; i_cut<n_cuts; ++i_cut) {

      // Read the column 'i_cut' that refers to the current line
      double value;
      sscanf(line, format[i_cut], &value);

      // If 'value' is not in the allowed range, we have to discard it.
      // To do so, we flag the line, and exit from the loop on the cuts.
      if( (value < cuts[2*i_cut]) || (value > cuts[2*i_cut+1]) ) {
        keep_line[i_line] = _FALSE_;
        
        // Debug
        // printf("The following line (%lld) was discarded:\n", i_line);
        // printf("    %s\n", line);
        
        break;
      }
    }

    // If none of the cuts above fired, then we have one more line to keep
    if (keep_line[i_line] == _TRUE_)
      ++(*n_lines_kept);

  }


  // =====================================
  // = Read the data from the input file =
  // =====================================
  // Rewind the file
  rewind(fp);
  
  // Now that we know the number of lines, we can allocate the memory
  *pc1 = (double*)malloc(sizeof(double)*(*n_lines_kept));
  *pc2 = (double*)malloc(sizeof(double)*(*n_lines_kept));
  *pc3 = (double*)malloc(sizeof(double)*(*n_lines_kept));
  *pc4 = (double*)malloc(sizeof(double)*(*n_lines_kept));
  *pc5 = (double*)malloc(sizeof(double)*(*n_lines_kept));
  *pc6 = (double*)malloc(sizeof(double)*(*n_lines_kept));
  
  long long skipped_lines=0;
  
  // Read data from the file
  for (i_line=0; i_line<*n_lines_read_from_file; i_line++) {

    // We need to store the information in columns 1 to 6.  The pointer to, e.g., column '3' is
    // *pc3.  The address of the element 'i' of column 3 is (*pc3)+i or, if you prefer the bracket
    // notation, (*pc3)[i].
    // 'lf' is the format specifier for a double precision number (read in 16 digits)
    char *fgets_output = fgets(line, LINE_MAX, fp);
    
    // Skip the lines that do not satisfy our cuts
    if (keep_line[i_line] == _FALSE_) {
      ++skipped_lines;
      continue;
    }
    
    // This line survives all of our cuts
    int line_idx = i_line - skipped_lines;
    sscanf(line,"%lf,%lf,%lf,%lf,%lf,%lf",
      (*pc1)+line_idx, (*pc2)+line_idx, (*pc3)+line_idx, (*pc4)+line_idx, (*pc5)+line_idx, (*pc6)+line_idx);

    // If we encounter the end of file, then it means that the file was shorter than the
    // 'n_lines_read_from_file' provided as an argument
    if (fgets_output==NULL) {
      *n_lines_read_from_file = i_line;
      break;
    }

  }  

  // Debug
  // long long i=0;
  // for (i=0; i<*n_lines_kept; i++) {
  //   printf("%lf,%lf,%lf,%lf,%lf,%lf\n", (*pc1)[i], (*pc2)[i], (*pc3)[i], (*pc4)[i], (*pc5)[i], (*pc6)[i]);
  // }
  
  // Close the file
  fclose(fp);
	
  return _SUCCESS_;
}



/// ============================
/// = read_six_column_csv_file =
/// ============================
int read_six_column_csv_file (
        char *file_name,
        long long *n_lines,
        // Will contain the data in arrays
        // Notation: pcn = pointer to column 'n'
        double **pc1, double **pc2, double **pc3, double **pc4, double **pc5, double **pc6
        ) {
  
  // Cycle variable
	long long i=0;
		
	// Open the data file
	FILE *fp = fopen(file_name,"r");
	if (fp == NULL) {
    ERRORPRINT("Couldn't open file %s.\n", file_name);
    return _FAILURE_;
  }
	  
  // Buffer that will contain a line of the file
  char line[LINE_MAX];

	// =========================================
	// = Count the number of lines in the file =
	// =========================================
  // If the user gave a number different from 0 as n_lines, then she wants only to 
  // consider the lines up to that number.  Otherwise, she wants all the lines in
  // the file, and we need to count them.
  if (*n_lines==0) {
  	while ( fgets(line, LINE_MAX, fp) != NULL )
  			(*n_lines)++;
  }

	// =====================================
	// = Read the data from the input file =
	// =====================================
	// Rewind the file
	rewind(fp);
	
	// Now that we know the number of lines, we can allocate the memory
	*pc1 = (double*)malloc(sizeof(double)*(*n_lines));
	*pc2 = (double*)malloc(sizeof(double)*(*n_lines));
	*pc3 = (double*)malloc(sizeof(double)*(*n_lines));
	*pc4 = (double*)malloc(sizeof(double)*(*n_lines));
	*pc5 = (double*)malloc(sizeof(double)*(*n_lines));
	*pc6 = (double*)malloc(sizeof(double)*(*n_lines));
  
  // Read data from the file
  for (i=0; i<*n_lines; i++) {

    // We need to store the information in columns 1 to 6.  The pointer to, e.g., column '3' is
    // *pc3.  The address of the element 'i' of column 3 is (*pc3)+i or, if you prefer the bracket
    // notation, (*pc3)[i]
    char *fgets_output = fgets(line, LINE_MAX, fp);
    sscanf(line,"%lf,%lf,%lf,%lf,%lf,%lf", (*pc1)+i, (*pc2)+i, (*pc3)+i, (*pc4)+i, (*pc5)+i, (*pc6)+i);
    
    // Apparently using fscanf is dangerous. See:
    // http://bytes.com/topic/c/answers/217715-fscanf-problem
    // http://www.gidnetwork.com/b-59.html
    // http://stackoverflow.com/questions/1252132/difference-between-scanf-and-fgets
    // fscanf(fp,"%lf,%lf,%lf,%lf,%lf,%lf\n", (*pc1)+i, (*pc2)+i, (*pc3)+i, (*pc4)+i, (*pc5)+i, (*pc6)+i);

    
    // If we encounter the end of file, then it means that the file was shorter than the
    // 'n_lines' provided as an argument
    if (fgets_output==NULL) {
      *n_lines = i;
      break;
    }

  }  
	
  // // Debug
  // for (i=0; i<*n_lines; i++) {
  //   printf("%lf,%lf,%lf,%lf,%lf,%lf\n", (*pc1)[i], (*pc2)[i], (*pc3)[i], (*pc4)[i], (*pc5)[i], (*pc6)[i]);
  // }

	// Close the file
	fclose(fp);
	
  return _SUCCESS_;
}


/// =====================================
/// = Function to read and split a file =
/// =====================================
int split_six_column_csv_file (char *in_filename, char *format, double side_length, int n_bins_grid) {
  
  // Cycle variable
  int i_x=0, i_y=0, i_z=0, i=0;
  
  // Counters
  long long n_lines=0, n_lines_in_blocks=0;
  
  // String for the output file name
  char output_file_name[NAME_MAX];
  
	// Temporary variables
  double x, y, z;
  
  // Array with the range of the grid
  double *grid_range = lin_space(0., side_length, n_bins_grid+1);
  double bin_width = side_length/((double)(n_bins_grid));
    
  // Avoid division by zero
  if (bin_width == 0.) {
    printf("ERROR, split_six_column_csv_file: stop to avoid division by zero for bin_width\n");
    return _FAILURE_;
  }

  // Some debug
  // printf("Divisions of the sides of the box:\n");
  // for (i=0; i<n_bins_grid+1; ++i)
  //   printf("%5f ", grid_range[i]);
  // printf("\n");

	// Open the input data file
	FILE *in_file = fopen(in_filename,"r");	
	if (in_file==NULL) {
    printf("ERROR, split_six_column_csv_file: could not open file %s.\n", in_filename);
    return _FAILURE_;
  }
	
	// ===============================
	// = Initialise the output files =
	// ===============================
  // Create the output directory
  mkdir (dirname(format), 0777);
	
	// Create an array of FILEs corresponding to the output files that will be generated
  int n_out_files = n_bins_grid*n_bins_grid*n_bins_grid;
  FILE **out_files = (FILE**)malloc(sizeof(FILE*)*n_out_files);
  
  // Create the output FILEs, and assign names to them in an automated way using
  // the 'format' parameter given by the user.
  for (i_x=0; i_x<n_bins_grid; ++i_x) {
    for (i_y=0; i_y<n_bins_grid; ++i_y) {
      for (i_z=0; i_z<n_bins_grid; ++i_z) {
        
        // Collapse the three grid indices to a single index that is used
        // to address the FILEs array
        int out_file_index = n_bins_grid*n_bins_grid*(i_x) + n_bins_grid*i_y + i_z;
        
        // The naming starts from bin 1-1-1, not from bin 0-0-0
        sprintf(output_file_name,format,i_x+1,i_y+1,i_z+1);
        out_files[out_file_index] = fopen(output_file_name,"w");
        
      }
    }
  }
  
  // =============
  // = Main loop =
  // =============
  // Buffer that will contain a line of the file
  char line[LINE_MAX];
  
  // Here we read in each line of the file and place it in the corresponding block of the grid
	while ( fgets(line, LINE_MAX, in_file) != NULL ) {
  
    // printf("Beginning of loop\n");      
    
    // Get the first three entries of the considered line
    sscanf(line,"%lf,%lf,%lf", &x, &y, &z);
    ++n_lines;
    
    // printf("Read file\n");
     
    // Here we use the fact that the grid is linearly spaced and that
    // we are binning it starting from x_min = y_min = z_min = 0.
    // It that was not the case, we should include minus x_min, y_min
    // and z_min factors (i.e. we should write x - x_min instead of just x,
    // and the same for y and z)
    i_x = (int)(x/bin_width);
    i_y = (int)(y/bin_width);
    i_z = (int)(z/bin_width);

    // printf("Done divisions.\n");  

    // If a galaxy is out of bounds, ignore it.
    if(  (i_x >= n_bins_grid) || (i_y >= n_bins_grid) || (i_z >= n_bins_grid)  ) 
      continue;

    // printf("Galaxy is in\n");
    // printf("i_x = %d,  i_y = %d,  i_z = %d\n", i_x, i_y, i_z);        
    
    // Write the line in one of the files if the line is within bounds
    int out_file_index = n_bins_grid*n_bins_grid*(i_x) + n_bins_grid*i_y + i_z;
    
    // printf("Writing to file %d\n", out_file_index);
    
    fputs(line, out_files[out_file_index]);
    n_lines_in_blocks+=1;
    
    // printf("End of loop\n");

  } // End while loop on file

  // Check that all the lines in the in_file were sorted out correctly
  if (n_lines!=n_lines_in_blocks)
    printf("WARNING: The input file had %lld lines while the output files have a total of %lld lines\n", n_lines, n_lines_in_blocks);

  // Close the output files 
  for (i=0; i<n_out_files; ++i)
    fclose(out_files[i]);
    
  // Close the input file
  fclose(in_file);
  free(out_files);
  
  return _SUCCESS_;
  
}



/// ===================
/// = Sorting indices =
/// ===================
// Definitions needed by the sorting functions
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp; 
#define M 7
#define NSTACK 50

int get_sorting_indices(double *array, long long n, long long *indx) {
  
  // Careful here with the indices. If your vectors are indexed from 1..n, change ir=n and l=1.
  long long i, indxt, ir=n-1, itemp, j, k, l=0; 
  long long jstack=0; 
  long long *istack;


  // Partitioning element in the array
  double a;

  // Allocate memory for the temporary array with NSTACK >= 2 log2 (n) elements. For the numbers we consider,
  // 50 is more than enough.
  istack=(long long*)malloc(sizeof(long long)*NSTACK);

  // Fill in the index array with integers from 0..n-1
  for (j=0;j<n;j++)
    indx[j]=j; 
  

  for (;;) {
  
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
        indxt=indx[j]; 
        a=array[indxt];
        for (i=j-1;i>=l;i--) {
          if (array[indx[i]] <= a) 
          break;
          indx[i+1]=indx[i];
        }
        indx[i+1]=indxt;
      }
      if (jstack == 0) 
      break; 
      ir=istack[jstack--]; 
      l=istack[jstack--];
    } 
    else {
      k=(l+ir) >> 1;
      SWAP(indx[k],indx[l+1])
      if (array[indx[l]] > array[indx[ir]]) {
        SWAP(indx[l],indx[ir])
      }
      if (array[indx[l+1]] > array[indx[ir]]) {
        SWAP(indx[l+1],indx[ir])
      }
      if (array[indx[l]] > array[indx[l+1]]) {
        SWAP(indx[l],indx[l+1])
      }
      i=l+1;
      j=ir; 
      indxt=indx[l+1];
      a=array[indxt]; 
    
      for (;;) {
        do i++; 
        while (array[indx[i]] < a); 
        do j--; 
        while (array[indx[j]] > a); 
        if (j < i) 
        break; 
        SWAP(indx[i],indx[j])
      }
      indx[l+1]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack > NSTACK) 
        printf("NSTACK too small in get_sorting_indices."); 
      if (ir-i+1 >= j-l) {
        istack[jstack]=ir; 
        istack[jstack-1]=i; 
        ir=j-1;
      } 
      else { 
        istack[jstack]=j-1; 
        istack[jstack-1]=l; 
        l=i;
      } 
    }
  }

  free(istack);

  return _SUCCESS_;
}

/// =================
/// = Sorting array =
/// =================
// Sort an array according to an array of integers.
int sort_array(double *array, long long n, long long *idx_array) {
  
  long long i=0;
  double *aux_array = (double*)malloc(sizeof(double)*n);
  
  for (i=0; i<n; ++i) {
    aux_array[i] = array[idx_array[i]];
  }
  
  for (i=0; i<n; ++i) {
    array[i]=aux_array[i];
  }
         
  free(aux_array);
     
  return _SUCCESS_;
  
}



/// ============================
/// = Calculate distances      =
/// ============================
// Given the coordinates of two 3D points, calculate the Euclidean distance between them.

double get_distance(double x1, double x2, double y1, double y2, double z1, double z2 ) {
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}




/// ===============================
/// = calculate infall velocities =
/// ===============================
// Given two galaxies, calculate the relative infall velocity.
double get_relative_infall_vel(
    double x1, double x2, double y1, double y2, double z1, double z2,
    double vx1, double vx2, double vy1, double vy2, double vz1, double vz2,
    double distance) {

  return ((x1-x2)*(vx1-vx2)+(y1-y2)*(vy1-vy2)+(z1-z2)*(vz1-vz2))/distance;

}



/// =============
/// = lin_space =
/// =============
// Return an array with linearly spaced points
double * lin_space ( double x_min, double x_max, int n_points) {

	int j;
  double *x = malloc(n_points*sizeof(double));

	x[0]=x_min;
  x[n_points-1]=x_max;
	double step = (x_max-x_min)/(n_points-1);

	for (j=1; j<n_points-1; ++j)
	    x[j] = x_min + j*step;
	    
	return x;
		
}


/// =============
/// = log_space =
/// =============
// Return an array with logarithmically spaced points.
double * log_space (double x_min, double x_max, int n_points) {

	int j;
  double *x = malloc(n_points*sizeof(double));

	x[0] = x_min;
  x[n_points-1]=x_max; 
	double step = pow (x_max/x_min, 1.0/(n_points-1));	

  for (j=1; j<n_points-1; ++j)  
    x[j] = x[j-1]*step;    

	return x;
	
}
