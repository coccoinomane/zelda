Zelda is a command-line tool to generate correlation functions and power spectra
from a galaxy catalogue.

Zelda has been tested on the data from the popular Millennium simulation 
(http://www.mpa-garching.mpg.de/galform/virgo/millennium/) and it led
to the publication of the following paper:
"Using galaxy pairs as cosmological tracers"
http://fr.arxiv.org/abs/1204.5761


DESCRIPTION

Zelda is written in C. Its structure is modular and flexible, and it was heavily
inspired by that of the cosmological Boltzmann code CLASS (http://class-code.net/).

Zelda is a parallel code, via the OpenMP standard. Therefore, it can use all the
cores in a single node, but it cannot span multiple nodes. Set the number of
cores you want to use on the command line via the environment variable
OMP_NUM_THREADS. For example, if you run Zelda on a laptop with 8 cores
on a bash shell, you might want to execute 'export OMP_NUM_THREADS=8' before
running Zelda.


INSTALLATION

Personalise the 'makefile' and run 'make all'.


SHORT USER GUIDE

The most important file for a new user is 'params_explanatory.ini'. It is a
text-only parameter file that can be fed to Zelda with the command
./zelda params_explanatory.ini
Read this documented file to learn about the most important parameters in Zelda.
The file can be also used as a template for creating your custom parameter files.

The directory structure of Zelda is important to learn how the code works:

* The 'source' directory contains the main source files in C. Each file
corresponds to a module in Zelda.

* The 'tool' directory contains accessory source files in C with
purely numerical functions or utility functions.

* The 'main' directory contains the main source files, i.e. the executable
files, including zelda.c.

* The 'python' directory contains Python scripts to launch Zelda, including
the Zelda wrapper (zelda.py), the batch script (zelda_script.py) and a
rebinning function.

* The 'include' directory contains the declaration files (.h) for all the
C files in the 'source', 'main' and 'tools' directories.

* The 'scripts' directory contains accessory script files in bash or 
other scripting languages. For example, to fetch catalogues from
remote servers (e.g. the Millennium simulation server) or handling
catalogue files.

* The 'test' directory contains executable programs to test the outputs
of Zelda.


CREDITS

We wish to thank Julien Lesgourgues, Thomas Tram and Diego Blas for creating
CLASS!


CONTACT

Please contact Guido W. Pettinari at guido.pettinari@gmail.com if you
need any help with the code.

