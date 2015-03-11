#ifndef _COMMON_
#define _COMMON_

/// ======================
/// = Standard libraries =
/// ======================
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libgen.h>           // dirname, basename
#include <sys/stat.h>         // stat, mkdir
// #include <sys/types.h>
// #include <unistd.h>

#ifdef _OPENMP
#include "omp.h"
#endif

/// ==================
/// = My definitions =
/// ==================
#define _SUCCESS_ 0
#define _FAILURE_ 1

#define _FALSE_ 0
#define _TRUE_ 1

#define _BADFLAG_ -999

/// =============================
/// = My preprocessor functions =
/// =============================
//// DEBUGGING MACROS
// Debugging macros so we can pin down message origin at a glance.  Note that
// the ... arguments cannot be empty, if there is nothing to pass use ""
#define STDERRPRINT(...)       fprintf(stderr, __VA_ARGS__)
#define WHERESTR  "[file %s, line %d]: "
#define WHEREARG  __FILE__, __LINE__
#define DEBUGPRINT(_fmt, ...)  STDERRPRINT(WHERESTR _fmt, WHEREARG, __VA_ARGS__)

//// ERROR MACROS
#define ERRORSTR  "ERROR, [%s, %d]: "
#define ERRORARG  __FILE__, __LINE__
#define ERRORPRINT(_fmt, ...)  STDERRPRINT(ERRORSTR _fmt, ERRORARG, __VA_ARGS__)

//// WARNING MACROS
#define WARNINGSTR  "WARNING, %s: "
#define WARNINGARG  __func__
#define WARNINGPRINT(_fmt, ...)  STDERRPRINT(WARNINGSTR _fmt, WARNINGARG, __VA_ARGS__)


/// ==================
/// = Error handling =
/// ==================
// All of these macros are taken from CLASS (http://lesgourg.web.cern.ch/lesgourg/class.html),
// so all the credit goes to them

#define _ERRORMSGSIZE_ 2048 /**< generic error messages are cut beyond this number of characters */
typedef char ErrorMsg[_ERRORMSGSIZE_]; /**< Generic error messages */

/* Function wrapper.  The do-while loop is needed because the preprocessor does not
   like trailing semicolons. */
#define zelda_call(function)						                                  \
  do {									                                                  \
    if (function == _FAILURE_) {					                                \
      return _FAILURE_;							                                      \
    }									                                                    \
  } while(0);

/* macro for testing condition and returning error if condition is true;
   args is a variable list of optional arguments, e.g.: args="x=%d",x 
   args cannot be empty, if there is nothing to pass use args="" */
#define zelda_test(condition, args...)				                      	    \
  do {									                                                  \
    if (condition) {                                                      \
      ErrorMsg Optional_arguments;				                              	\
      sprintf(Optional_arguments,args);			                          		\
      ERRORPRINT("Condition (%s) is true; %s\n",                          \
        #condition, Optional_arguments);           			                	\
      return _FAILURE_;							                                      \
    }									                                                    \
  } while(0);


/* macro for opening file and returning error if it failed */
#define zelda_open(pointer,						                                    \
		   filename,						                                              \
  	   mode)        				                                              \
  do {									                                                  \
    pointer=fopen(filename,mode);				                                	\
    if (pointer == NULL) {					                                    	\
      ERRORPRINT("could not open %s with name %s and mode %s",	          \
	      #pointer,filename,#mode);				                                	\
      return _FAILURE_;						                                      	\
    }									                                                    \
  } while(0);

/* macro for allocating memory and returning error if it failed */
#define zelda_alloc(pointer, size)                                        \
  do {									                                                  \
    pointer=malloc(size);						                                      \
    if (pointer == NULL) {						                                    \
      int size_int;							                                          \
      size_int=size;							                                        \
      ERRORPRINT("Could not allocate %s with size %d",		                \
	      #pointer,size_int);					                                      \
      return _FAILURE_;							                                      \
    }									                                                    \
  } while(0);

/* macro for allocating memory and returning error if it failed */
#define zelda_calloc(pointer, num, size)                                  \
  do {									                                                  \
    pointer=calloc(num,size);				                                      \
    if (pointer == NULL) {						                                    \
      int size_int;							                                          \
      size_int=size;							                                        \
      ERRORPRINT("Could not allocate %s with size %d",		                \
	      #pointer,size_int);					                                      \
      return _FAILURE_;							                                      \
    }									                                                    \
  } while(0);

// ============================
// = Enums & type definitions =
// ============================
// PATH_MAX = 1024 on a Macbook Pro with OSX 10.7 (Lion)
typedef char FileName[PATH_MAX];




#endif