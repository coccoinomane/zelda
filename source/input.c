/** @file input.c Documented input module.
 *
 * Julien Lesgourgues, 27.08.2010    
 */

#include "input.h" 

/**
 * Use this routine to extract initial parameters from files 'xxx.ini'
 * and/or 'xxx.pre'. They can be the arguments of the main() routine.
 *
 * If class is embedded into another code, you will probably prefer to
 * call directly input_init() in order to pass input parameters
 * through a 'file_content' structure.
 */

int input_init_from_arguments(
			      int argc, 
			      char **argv,
       		  params_struct * params,
       		  data_struct * data,
  		  		ErrorMsg errmsg
			      ) {

  /** Summary: */

  /** - define local variables */

  struct file_content fc;
  struct file_content fc_input;
  struct file_content fc_precision;

  char input_file[_ARGUMENT_LENGTH_MAX_];
  char precision_file[_ARGUMENT_LENGTH_MAX_];

  int i;
  char extension[5];

  /** - Initialize the two file_content structures (for input
      parameters and precision parameters) to some null content. If no
      arguments are passed, they will remain null and inform
      init_params() that all parameters take default values. */

  fc.size = 0;
  fc_input.size = 0;
  fc_precision.size = 0;
  input_file[0]='\0';
  precision_file[0]='\0';

  /** If some arguments are passed, identify eventually some 'xxx.ini'
      and 'xxx.pre' files, and store their name. */

  if (argc > 1) {
    for (i=1; i<argc; i++) {
      strncpy(extension,(argv[i]+strlen(argv[i])-4),4);
      extension[4]='\0';
      if (strcmp(extension,".ini") == 0) {
	zelda_test(input_file[0] != '\0',
		   "You have passed more than one input file with extension '.ini', choose one.");
	strcpy(input_file,argv[i]);
      }
      if (strcmp(extension,".pre") == 0) {
	zelda_test(precision_file[0] != '\0',
		   "You have passed more than one precision with extension '.pre', choose one.");
	strcpy(precision_file,argv[i]);
      }
    }
  }
  
  /** - if there is an 'xxx.ini' file, read it and store its content. */

  if (input_file[0] != '\0')
    
    zelda_call(parser_read_file(input_file,&fc_input,errmsg));
    
  /** - if there is an 'xxx.pre' file, read it and store its content. */

  if (precision_file[0] != '\0')
    
    zelda_call(parser_read_file(precision_file,&fc_precision,errmsg));
    
  /** - if one or two files were read, merge their contents in a
      single 'file_content' structure. */

  if ((input_file[0]!='\0') || (precision_file[0]!='\0'))
    zelda_call(parser_cat(&fc_input,&fc_precision,&fc,errmsg));

  zelda_call(parser_free(&fc_input));
  zelda_call(parser_free(&fc_precision));
  
  /** - now, initialize all parameters given the input 'file_content'
      structure.  If its size is null, all parameters take their
      default values. */

  zelda_call(input_init(
      &fc,
			params,
      data,
      errmsg));

  zelda_call(parser_free(&fc));

  return _SUCCESS_;
}

/**
 * Initialize each parameters, first to its default values, and then
 * from what can be interpreted from the values passed in the input
 * 'file_content' structure. If its size is null, all parameters keep
 * their default values.
 */
// 
int input_init(
        struct file_content * pfc,  
        params_struct * params,
        data_struct * data,
        ErrorMsg errmsg		 
        ) {

  /** Summary: */

  /** - define local variables */

  int flag1,flag2,flag3;
  double param1,param2,param3;
  int n,entries_read;
  int int1,fileentries;
  double fnu_factor;
  double * pointer1;
  char string1[_ARGUMENT_LENGTH_MAX_];

  double Omega_tot;

  int i;

  FILE * param_output;
  FILE * param_unused;
  char param_output_name[_LINE_LENGTH_MAX_];
  char param_unused_name[_LINE_LENGTH_MAX_];


  /** - set all parameters (input and precision) to default values */

  zelda_call(input_default_params(
        params,
        data));



  /** - if entries passed in file_content structure, carefully read
      and interpret each of them, and tune accordingly the relevant
      input parameters */

  // ============
  // = Modality =
  // ============
  zelda_call(parser_read_string(pfc,
       "modality",
       &(string1),
       &(flag1),
       errmsg));

  if (flag1 == _TRUE_) {
    int flag2 = _FALSE_;
    if (strcmp(string1,"vinfall") == 0) {
      params->modality = VINFALL;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"vinfall_isolation") == 0) {
      params->modality = VINFALL_ISOLATION;
      flag2=_TRUE_;
    }
    zelda_test(flag2 == _FALSE_,
      "could not identify what to compute, set modality to a value between the following: vinfall, vinfall_isolation");
  }
  
  
  // ==========================
  // = File naming parameters =
  // ==========================
  zelda_read_string("root",params->root);
  zelda_read_string("input_format", params->input_format);
  zelda_read_string("results_filename", params->results_filename);

  // ===========================
  // = Cuts-related parameters =
  // ===========================
  parser_read_list_of_integers(
        pfc,
        "columns_to_cut",
        &entries_read,
        &(params->columns_to_cut),
        &flag2,
        errmsg
        );
  if (flag2==_TRUE_)
    params->n_cuts = entries_read;

  parser_read_list_of_doubles(
        pfc,
				"cuts",
				&entries_read,
				&(params->cuts),
				&flag2,
				errmsg);
  if (params->n_cuts>0)								
    zelda_test((entries_read!=2*params->n_cuts) || flag2 == _FALSE_, "the 'cuts' field should have exactly double the entries of the field 'columns_to_cut'");
				
  // ======================
  // = Separation binning =
  // ======================
  zelda_read_double("r_min", params->r_min);
  zelda_read_double("r_max", params->r_max);
  zelda_read_int("n_box_bins", params->n_box_bins);  

  // Parse the binning mode
  zelda_call(parser_read_string(pfc,
       "binning_mode",
       &(string1),
       &(flag1),
       errmsg));

  if (flag1 == _TRUE_) {
    int flag2 = _FALSE_;
    if (strcmp(string1,"lin") == 0) {
      params->binning_mode = LIN_BINNING;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"log") == 0) {
      params->binning_mode = LOG_BINNING;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"custom") == 0) {
      params->binning_mode = CUSTOM_BINNING;
      flag2=_TRUE_;
    }    
    zelda_test(flag2 == _FALSE_,
      "could not identify the requested separation binning, set 'binning_mode' to a value between the following: lin, log, custom");
  }
  

  // If the user did not provide a custom binning, then we read the number
  // of separation bins from the parameters file.
  if (params->binning_mode != CUSTOM_BINNING) {
    zelda_read_int("n_separation_bins",params->n_separation_bins);
  }
  // It the user selected to provide a custom binning, read it from the 'custom_bin_edges'
  // and 'custom_bin_midpoints' fields.
  else {

    // We read the edges of the bins and 'n_separation_bins' from the field 'custom_bin_edges'
    parser_read_list_of_doubles(
          pfc,
  				"custom_bin_edges",
  				&entries_read,
  				&(params->r_edges),
  				&flag2,
  				errmsg);

  	// Check that the user provided at least two bin edges
    zelda_test((entries_read<2) || (flag2 == _FALSE_), "the 'custom_bin_edges' field should have at least two entries.");

    // Define n_separation_bins, r_min, r_max
    params->n_separation_bins = entries_read-1;
    params->r_min = params->r_edges[0];
    params->r_max = params->r_edges[params->n_separation_bins];    

    // Check that the bins are sorted in increasing order
    int i_bin;
    for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin)
      zelda_test(params->r_edges[i_bin] > params->r_edges[i_bin+1],
        "please ensure that the field 'custom_bin_edges' is a comma separated list sorted in increasing order"); 

    // Read the midpoints
    parser_read_list_of_doubles(
          pfc,
  				"custom_bin_midpoints",
  				&entries_read,
  				&(params->r_midpoints),
  				&flag2,
  				errmsg);

    // If no midpoints are specified, compute them linearly from the r_edges array.
    if (flag2 == _FALSE_) {
      free(params->r_midpoints);
      zelda_alloc(params->r_midpoints, params->n_separation_bins*sizeof(double));
      for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin)
        params->r_midpoints[i_bin] = 0.5*(params->r_edges[i_bin+1] + params->r_edges[i_bin]);
    }
    // If the midpoints are specified, check that they are actually midpoints.
    else if (entries_read==params->n_separation_bins) {  
      for (i_bin=0; i_bin<params->n_separation_bins; ++i_bin) {
        double r = params->r_midpoints[i_bin];
        zelda_test(  (r >= params->r_edges[i_bin+1]) || (r <= params->r_edges[i_bin]),
          "The midpoints you specified in 'custom_bin_midpoints' are not in the correct range wrt 'custom_bin_edges'","");
      } 
    }
    else {
      ERRORPRINT("Please make sure that 'custom_bin_midpoints' is a comma separated list with one element less than 'custom_bin_edges'.\n","");
      return _FAILURE_;
    }
      
  }  // End of if (CUSTOM_BINNING)




	// ======================
	// = Intermediate files =
	// ======================
  zelda_call(parser_read_string(pfc,
       "store_blocks_results",
       &(string1),
       &(flag1),
       errmsg));       
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL) || atoi(string1) > 0))
    params->store_blocks_results = 1;

	// Output format for blocks' results
  zelda_read_string("output_format", params->output_format);

  // =======================
  // = Isolation criterion =
  // =======================
  // If r_isolation is not specified in the parameters file,
  // just take it to be equal to r_max. THIS MUST BE BELOW zelda_read_double(R_MAX)!
  zelda_call(parser_read_double(pfc,"r_isolation",&param1,&flag1,errmsg));
  if (flag1 == _FALSE_)
    params->r_isolation = params->r_max;
  else
    params->r_isolation = param1;

  zelda_read_int("n_neighbours_max",params->n_neighbours_max);  
  
  // ================
  // = Edge effects =
  // ================
  zelda_read_double("side_of_box", params->side_of_box);
  zelda_call(parser_read_string(pfc,
       "remove_galaxies_on_the_edges",
       &(string1),
       &(flag1),
       errmsg));       
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL) || atoi(string1) > 0))
    params->remove_galaxies_on_the_edges = _TRUE_;
  else
    params->remove_galaxies_on_the_edges = _FALSE_;
  
  // ====================
  // = Debug parameters =
  // ====================
  zelda_read_int("n_gals_per_block",params->n_gals_per_block);  
  zelda_read_int("verbose",params->verbose);    


  // =============================
  // = Dump parameters to a file =
  // =============================
  zelda_call(parser_read_string(pfc,"write parameters",&string1,&flag1,errmsg)); 

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {

    sprintf(param_output_name,"%s%s",params->root,"parameters.ini");
    sprintf(param_unused_name,"%s%s",params->root,"unused_parameters");

    zelda_open(param_output,param_output_name,"w");
    zelda_open(param_unused,param_unused_name,"w");

    fprintf(param_output,"# List of input/precision parameters actually read\n");
    fprintf(param_output,"# (all other parameters set to default values)\n");
    fprintf(param_output,"#\n");
    fprintf(param_output,"# This file, written by Zelda, can be used as the input file\n");
    fprintf(param_output,"# of another run\n");
    fprintf(param_output,"#\n");

    fprintf(param_unused,"# List of input/precision parameters passed\n");
    fprintf(param_unused,"# but not used (just for info)\n");
    fprintf(param_unused,"#\n");

    for (i=0; i<pfc->size; i++) {
      if (pfc->read[i] == _TRUE_)
        fprintf(param_output,"%s = %s\n",pfc->name[i],pfc->value[i]);
      else
        fprintf(param_unused,"%s = %s\n",pfc->name[i],pfc->value[i]);
    }
    fprintf(param_output,"#\n");

    fclose(param_output);
    fclose(param_unused);
  }

  return _SUCCESS_;

}

/** 
 * All default parameter values (for input parameters)
 */ 

int input_default_params(
      params_struct * params,
 		  data_struct * data
      ) {

  params->modality = VINFALL;
  strcpy(params->root, "./");
  strcpy(params->input_format, "gal_blocks_%d_%d_%d.dat");
  strcpy(params->results_filename, "results_font2008_z2.dat");
  params->n_box_bins=10;
  params->side_of_box=500.;
  params->remove_galaxies_on_the_edges=_TRUE_;
  params->n_gals_per_block=0;
	params->store_blocks_results=_FALSE_;
  strcpy(params->output_format, "vinfall_blocks_%d_%d_%d.dat");
  params->verbose=0;
  params->n_neighbours_max=1;
    
  // Binning related parameters
  params->r_min=0.;
  params->r_max=50.;   // r_max also determines the default value of r_isolation
  params->n_separation_bins=10;
  params->binning_mode=LIN_BINNING;
  
  // Cuts-related parameters
  params->n_cuts=0;
  params->columns_to_cut=NULL;
  params->cuts=NULL;

  return _SUCCESS_;

}

