## Makefile for the Zelda project
## Created by Guido W. Pettinari on 21/10/2011 based
## on the makefile for CLASS (http://class-code.net/).
## Last modified by Guido W. Pettinari on 31/10/2011

CC       						= gcc
LDFLAGS  						= -lgomp
INCLUDES 						= include
OPTIMIZATION_FLAGS 	= -O3 -ffast-math
OPENMP   						= -fopenmp
CFLAGS   						= -g $(OPTIMIZATION_FLAGS) $(OPENMP)

## Rule to create the building directory (i.e. the directory where
## all the temporary .o files will be stored). The file .base is just
## a placeholder that prevents make to create the build directory 
## everytime.
WRKDIR = build
.base:
	if ! [ -a $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base
vpath .base build

## Tell make where to look for files
vpath %.c source:tools:main:tests
vpath %.o build

## Rule for compiling a .c file
%.o:  %.c .base
	$(CC) $(CFLAGS) -I$(INCLUDES) -c $< -o $(WRKDIR)/$*.o

## Object files
TOOLS = tools.o parser.o
INPUT = input.o
PARAMS = params.o
DATA_TOOLS = data_tools.o data_block_tools.o

## Mains
ZELDA = zelda.o
FILE_SPLITTER = file_splitter.o
PRINT_PARAMS = print_params.o
TEST_CUTS = test_cuts.o

## Rules
all: zelda file_splitter print_params

zelda: $(DATA_TOOLS) $(PARAMS) $(TOOLS) $(INPUT) $(ZELDA)
	$(CC) $(LDFLAGS) -o $@ $(addprefix $(WRKDIR)/,$(notdir $^))

file_splitter: $(TOOLS) $(FILE_SPLITTER)
	$(CC) $(LDFLAGS) -o $@ $(addprefix $(WRKDIR)/,$(notdir $^))

print_params: $(DATA_TOOLS) $(PARAMS) $(TOOLS) $(INPUT) $(PRINT_PARAMS)
	$(CC) $(LDFLAGS) -o $@ $(addprefix $(WRKDIR)/,$(notdir $^))

test_cuts: tools.o $(TEST_CUTS)
	$(CC) $(LDFLAGS) -o $@ $(addprefix $(WRKDIR)/,$(notdir $^))

clean: .base
	rm -rf $(WRKDIR);