#include "zelda.h"


int main (int argc, char ** argv) {

  // =========================
  // = Create the structures =
  // =========================
  params_struct params;
  data_struct data;
  ErrorMsg error_message;


  // =======================
  // = Fill the structures =
  // =======================
  if( input_init_from_arguments(argc, argv, &params, &data, error_message) == _FAILURE_ ) {
    printf("Something went wrong when loading 'input_init_from_arguments'.\n");
    return _FAILURE_;
  }

  if( params_init(&params) == _FAILURE_ ) {
    printf("Something went wrong when loading 'params_init'.\n");
    return _FAILURE_;
  }

  if( data_init(&params, &data) == _FAILURE_ ) {
    printf("Something went wrong when loading 'data_init'.\n");
    return _FAILURE_;
  }
  

  // =======================
  // = Free the structures =
  // =======================
  if (params_free (&params) == _FAILURE_) {
    printf("Something went wrong when loading 'params_free'.\n");    
    return _FAILURE_;  
  }

  if (data_free (&data) == _FAILURE_) {
    printf("Something went wrong when loading 'data_free'.\n");        
    return _FAILURE_;  
  }
    

  
  
  return _SUCCESS_;

}



