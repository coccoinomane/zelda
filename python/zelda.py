#!/usr/bin/env python

"""zelda.py: The Zelda package is a Python wrapper for the galaxy pair C code Zelda."""

__author__      = "Guido Walter Pettinari"
__copyright__   = "Copyright 2012, Planet Earth"

import os
import sys
import re
from shutil import copy
import fileinput
import csv
import numpy as np
import shutil
import fileinput
from time import gmtime, strftime
import re
import inspect
def func_name():
    return inspect.stack()[1][3]

import errno
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise


class zelda_run:
    
    """Dictionary that contains the parameters that will be used to run Zelda"""
    params = None
    """File name of the parameter file that will be created"""
    params_filename = None
    """String that identifies this run"""
    identifier = None
    """Standard template to construct a parameters file"""
    params_template =  {
                       "modality"							:	"vinfall",
                       "columns_to_cut"						:	"",
                       "cuts"								:	"",
                       "input_format"						:	"",
                       "output_format"						:	"./vinfall_blocks/vinfall_blocks_%d_%d_%d.dat",
                       "n_box_bins"							: 	1,
                       "results_filename"					: 	"",
                       "r_isolation"						: 	"",
                       "n_neighbours_max"                   :   1,
                       "r_min"								: 	1e-3,
                       "r_max"								: 	70,
                       "n_separation_bins"					:	50,
                       "binning_mode"    					: 	"lin",
                       "n_gals_per_block"					: 	0,
                       "side_of_box"						: 	500,
                       "remove_galaxies_on_the_edges"		:	"y",
					   "store_blocks_results"				: 	"n",
                       "verbose"							:	3
                       }

    """Additional text to indentify the run, that will be appended to self.identifier.
    It should be a string.  The line below is just an example."""
    extra_naming = "guo2010a_z0"
                                    

    def __init__(self, zelda_bin, out_dir="./", custom_params=None, extra_naming=None):
        """Initialise a Zelda object with some custom parameters specified via a dictionary"""
        self.zelda_bin = zelda_bin
        self.out_dir = out_dir
        self.params = self.params_template
        if isinstance(custom_params, dict):
            self.update(custom_params)
        if isinstance(extra_naming, basestring):
            self.extra_naming=extra_naming

    def update(self, custom_params, extra_naming=None):
        self.params.update(custom_params)
        if isinstance(extra_naming, basestring):
            self.extra_naming = extra_naming
        self.build_identifier()

    def build_identifier(self):
        """Construct a string that identifies this run and can be used to name
        output files accordingly"""
        # Example:
        # results_guo2010_z0_log_r1e-3to4_nbins15_total_isolation_single_block_no_galsonedges.dat
     
        # Build the identifier
        identifier_parts = [
            "method_part"       ,   self.params["modality"],
            "extranaming_part"  ,   self.extra_naming,
            "linlogcustom_part" ,   self.params["binning_mode"],
            "range_part"        ,   "r" + str(self.params["r_min"]) + "to" + str(self.params["r_max"]),
            "nbins_part"        ,   "nbins" + str(self.params["n_separation_bins"]),
            "ngrid_part"        ,   "ngrid" + str(self.params["n_box_bins"])
        ]

        if (self.params["modality"]=="vinfall_isolation"):
            identifier_parts.append("neighbours_part")
            identifier_parts.append("neigh" + str(self.params["n_neighbours_max"]))
            
        self.identifier = "_".join(identifier_parts[1::2])

        # Strip double underscores that may be there if any of the filename parts was empty
        self.identifier = re.sub("__", "_", self.identifier)
        
        # Identify the run as a test run if it does not analyse of the full dataset
        self.is_test_run = (self.params["n_gals_per_block"]!=0)
        
        # ===================
        # = Build filenames =
        # ===================
        # The root directory where everything will be stored
        ## Uncomment if you want a sub directory for every run
        # self.root = os.path.join(self.out_dir, self.identifier)
        ## Uncomment if you do not want a directory stucture
        self.root = self.out_dir        
        # Build the results filename according to the identifier
        self.params["results_filename"] = os.path.join(self.root, "results", self.identifier + ".dat")
        self.results_filename = self.params["results_filename"]        
        # Build the params filename according to the indentifier
        self.params_filename = os.path.join(self.root, "params", self.identifier + ".ini")
        # Build the log filename according to the indentifier
        self.log_filename = os.path.join(self.root, "logs", self.identifier + ".log")
    
    def dump_params_to_file(self):
        """Create a parameter file readable by Zelda using the values in self.params.
        Also create the directory structure for the results."""

        # Build the string that will identify this run, and create the directories
        if self.identifier == None:
            self.build_identifier()
        mkdir_p(os.path.dirname(self.params["results_filename"]))
        mkdir_p(os.path.dirname(self.params_filename))
        mkdir_p(os.path.dirname(self.log_filename))
                
        # Write to file the parameters
        fp = open(self.params_filename, 'w')
        fp.write("## File generated by %s on %s\n" % (sys.argv[0], strftime("%Y-%m-%d %H:%M:%S", gmtime())))
        for (k,v) in self.params.items():
            fp.write(k + " = " + str(v) + '\n')
        fp.close()
        print func_name() + ": generated parameter file '%s'" % self.params_filename
                
    
    def run(self):
        """Run Zelda using the parameters contained in self.params"""

        # First of all, we need to generate the string identifier of the run.
        # The filenames for the parameters, for the results, and for the log
        # are also created here.
        self.build_identifier()

        # ===========================================
        # = Check whether the result already exists =
        # ===========================================
        if (os.path.isfile(self.log_filename) & os.path.isfile(self.results_filename)):

            # Look in the log file for the string "Execution finished"
            with open(self.log_filename) as log_file:
                for line in log_file:
                    if "Execution finished" in line:
                        print "*** zelda_run.run: NOT RUNNING BECAUSE RESULT FILE ALREADY EXISTS IN '%s'" % self.results_filename
                        return 0

        # Create the parameter file to feed Zelda with.
        self.dump_params_to_file()
            
        # =============
        # = Run Zelda =
        # =============
        # Create the log file, and write the current time on it
        with open(self.log_filename,'w') as log_file:
            self.time_start = strftime("%Y-%m-%d %H:%M:%S", gmtime())            
            log_file.write("*** Execution started  on %s\n" % self.time_start)
        
        # Build the Zelda command and execute it. All output is redirected to the log file
        cmd = self.zelda_bin + " " + self.params_filename + " >>  " + self.log_filename + " 2>&1"
        print "Executing command '" + cmd + "'..."
        return_code = os.system(cmd)
        self.is_successful_run = (return_code==0)
        
        # Print a warning is something went wrong
        if (not(self.is_successful_run)):
            print "WARNING: Zelda exited with error code %d" % return_code
            with open(self.log_filename,'a') as log_file:
                self.time_end = strftime("%Y-%m-%d %H:%M:%S", gmtime())            
                log_file.write("*** WARNING: Zelda exited with code %d on %s\n" % (return_code, self.time_end))
                return 1

        # If the run was a success, append information on the execution time to the log_file.
        # This is done only if the run was succesful and if it was not a test run.
        if self.is_successful_run & (not(self.is_test_run)):
            print "The Zelda run has exited with SUCCESS."
            with open(self.log_filename,'a') as log_file:
                self.time_end = strftime("%Y-%m-%d %H:%M:%S", gmtime())            
                log_file.write("*** Execution finished on %s\n" % self.time_end)
                return 0

class zelda_data:

    """String that encodes how to read the blocks associated to a catalog"""
    input_format = None
    """List of integers that contains the columns that we are using to select the data"""
    columns_to_cut = None
    """List of floats with the minimum and maximum of each column to cut"""
    cuts = None

    def __init__(self):
        pass


def generate_params(input_format, isolation=False, r_iso=None, n_neighbours_max=None, r_min=None,
           r_max=None, verbose=None, binning_mode=None, n_separation_bins=None, n_grid=None,
           output_format=None, n_gals_per_block=None
           ):

    params_default = {
                     "modality"							:	"vinfall",
                     "columns_to_cut"					:	"",
                     "cuts"								:	"",
                     "input_format"						:	"",
                     "output_format"					:	"./vinfall_blocks/vinfall_blocks_%d_%d_%d.dat",
                     "n_box_bins"						: 	1,
                     "results_filename"					: 	"",
                     "r_isolation"						: 	"",
                     "n_neighbours_max"                 :   1,                     
                     "r_min"							: 	1e-3,
                     "r_max"							: 	70,
                     "n_separation_bins"				:	50,
                     "binning_mode"  					: 	"lin",
                     "n_gals_per_block"					: 	0,
                     "side_of_box"						: 	500,
                     "remove_galaxies_on_the_edges"		:	"y",
					 "store_blocks_results"				: 	"n",
                     "verbose"							:	3
                     }

    params = params_default

    if isinstance(input_format, basestring):
        params["input_format"] = input_format
    if isinstance(n_grid, int):
        params["n_box_bins"] = n_grid
    if isinstance(binning_mode, basestring):
        params["binning_mode"] = binning_mode
    if isinstance(n_separation_bins, int):
        params["n_separation_bins"] = n_separation_bins
    if isinstance(r_max, (int,float)):
        params["r_max"] = r_max
    if isinstance(n_gals_per_block, int):
        params["n_gals_per_block"] = n_gals_per_block
    if isinstance(verbose, int):
        params["verbose"] = verbose
    
    if isolation:
        params.update({
            "modality"          :   "vinfall_isolation",
            "r_max"             :   r_iso,
            "n_neighbours_max"  :   n_neighbours_max,
            "n_box_bins"        :   n_grid,
            "remove_galaxies_on_the_edges": "y"
            })
                
    return params