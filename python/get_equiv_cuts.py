#!/usr/bin/env python

import os
import sys
import re
from shutil import copy
import fileinput
import csv
import numpy as np
# import pylab
import fileinput
from time import gmtime, strftime
import re


# Function that gets a reference column and finds the number of elements that fullfil the
# cut (N). Then, it sorts the other column to cut, finds the first (or last, according to the sorting)
# N elements and gets the cut for that column from them. The mode argument is used to set the proportionality
# of the columns. If mode=0, the columns compared are proportional while if mode=1, they are inversely so.
def get_column_cut(ref_column, cut_ref_column, column_to_cut, mode):
  # Get the reference column into an array
  # ref_col = np.loadtxt(filename, delimiter=',', usecols=[ref_column-1])
  # Calculate the number of survivors
  number_survivors = np.sum(ref_column >= cut_ref_column)
  # Load the column to get the cut from
  # col_to_cut= np.loadtxt(filename, delimiter=',', usecols=[column_to_cut-1])
  # Sort the column to be cut (we introduce the absolute value so that we don't run into problems with 
  # negative values)
  sorted_col = np.sort(column_to_cut)
  if mode==0:
    # Reverse the sorted column
    sorted_col = sorted_col[::-1]
  # Now we select the value of the column corresponding to the cut we made
  return (sorted_col[number_survivors-1], number_survivors)



filename = sys.argv[1]

# Columns for np, stellarMass, rMag respectively
cols = np.array([7,9,10])

# Allocate array to contain data
max_nlines = 30000000
data = np.ndarray([max_nlines,3])

# Read in data
n_lines = 0
csv_reader = csv.reader(open(filename, "r"))
for line in csv_reader:
    data[n_lines,0] = np.float32(line[7])
    data[n_lines,1] = np.float32(line[9])
    data[n_lines,2] = np.float32(line[10])
    n_lines = n_lines+1    
    # Alas, the 'map' approach is much slower
    # data[n_lines] = np.array(map(np.float32,line))[cols]

print "**** Number of galaxy in the file"
print n_lines


# The np.loadtxt approach will suck up all the memory of the computer
# data = np.loadtxt(test_filename, delimiter=',', usecols=cols)


# The number of particle is always larger than 20 in the Millennium simulation
cuts_np = [100,500,1000,5000]
cuts_stellar_mass=[]
cuts_rmag=[]
n_survivors=[]
fraction_survivors=[]

# Assign names to the columns
np_col, stellarmass_col, rmag_col = [0,1,2]

for cut_np in cuts_np:
    
    # Find the equivalent cut for stellarmass and append it to the corresponding list
    (equivalent_cut, n_surviv) = get_column_cut(
                                            data[:n_lines,np_col],
                                            cut_np,
                                            data[:n_lines,stellarmass_col],
                                            0 # do not reverse order
                                            )

    cuts_stellar_mass.append(equivalent_cut)

    # Find the equivalent cut for rmag and append it to the corresponding list
    (equivalent_cut, survivors_number) = get_column_cut(
                                            data[:n_lines,np_col],
                                            cut_np,
                                            data[:n_lines,rmag_col],
                                            1 # reverse order
                                            )
    cuts_rmag.append(equivalent_cut)
    
    # The number of survivors is the same for both cuts
    n_survivors.append(n_surviv);

# We don't need data anymore
del data

print "**** Number of survivors"
print np.array(n_survivors)

print "**** Fraction of survivors"
print np.asarray(n_survivors,np.float64)/np.float64(n_lines)

print "**** Original cut in np"
print cuts_np

print "**** Equivalent cuts in stellarMass"
print cuts_stellar_mass

print "**** Equivalent cuts in rMag"
print cuts_rmag

