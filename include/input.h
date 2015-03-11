/** @file input.h Documented includes for input module */

#ifndef __INPUT__
#define __INPUT__

#include "common.h"
#include "parser.h"
#include "params.h"
#include "data_tools.h"
#include "data_block_tools.h"

/* macro for reading parameter values with routines from the parser */
#define zelda_read_double(name,destination)			                     	\
  do {						                                              			\
    zelda_call(parser_read_double(pfc,name,&param1,&flag1,errmsg));   \
    if (flag1 == _TRUE_)						                                  \
      destination = param1;						                                \
  } while(0);


#define zelda_read_int(name,destination)			                      	\
  do {									                                              \
    zelda_call(parser_read_int(pfc,name,&int1,&flag1,errmsg));        \
    if (flag1 == _TRUE_)					                                  	\
      destination = int1;				                                    	\
  } while(0);

#define zelda_read_string(name,destination)			                    	\
  do {								                                                \
    zelda_call(parser_read_string(pfc,name,&string1,&flag1,errmsg));  \
    if (flag1 == _TRUE_)						                                  \
      strcpy(destination,string1);					                          \
  } while(0);




/**************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int input_init_from_arguments(
		 int argc, 
		 char **argv,
		 params_struct * params,
		 data_struct * data,
 		 ErrorMsg errmsg
		 );

  int input_init(
		 struct file_content * pfc,  
		 params_struct * params,
		 data_struct * data,
		 ErrorMsg errmsg		 
		 );

  int input_default_params(
		 params_struct * params,
		 data_struct * data
	   );


#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
