#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
import re
import os

def flatten( list ):
    """Turn a list of lists into a list with the elements of the sublists.
    Trick taken from http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
    """
    return sum(list,[])


def rebin( xx, yy, nn, n_min, var, edges=None):
    """Rebin an array imposing a minimum number of points in each bin.
    
    Input:
    xx, yy, nn -> 1D array containing, respectively, the bin-midpoints, the values of the 
        bins and the number of data points for each bin.
           
    n_min -> minimum number of data points each bin must contain in the returned arrays.
    
    var -> 1D array with the variance associated to each bin.
    
    edges -> 1D array with the edges associated to the original binning.  If present, the
        function will return also the re-binned edges.
    
    Output:
    (xx, yy, nn, variance, sigma, err, [edges])
    
    """
    
    # Check that we have enough data points so that can make at least one bin
    assert np.sum(nn) > n_min
    n_total = np.sum(nn)    
    
    # Check that the input is consistent
    assert len(xx) == len(yy) == len(nn)
    n_bins = len(xx)
    
    # Understand whether the user wants to compute the edges
    edges_requested = isinstance(edges, list) or isinstance(edges, np.ndarray)
    
    # Work with arrays
    xx = np.asarray(xx, dtype=np.float64).copy()
    yy = np.asarray(yy, dtype=np.float64).copy()
    nn = np.asarray(nn, dtype=np.float64).copy()
    var = np.asarray(var, dtype=np.float64).copy()
    if edges_requested:
        edges = np.asarray(edges, dtype=np.float64).copy()
    
    # Cycle through the underpopulated bins, and find out with which bins they
    # need to be merged with to obtain a bin with more than 'n_min' data points.
    # All underpopulated bin will be merged with the next bin at their right. If this is
    # not enough to obtain 'n_min' data points, another merger with the next bin
    # will take place.  The process is iterated until the n>n_min requirement is met.
    bad_bins_idxs = (nn < n_min).nonzero()
    bins_to_merge = [];
    bins_to_delete = []
    for i in bad_bins_idxs[0]:
        
        # The 'bins_to_merge' variable keeps track of the bins that need to be merged.
        # The bin we are considering will already be in the list if the previous
        # bin was an underpopulated one. If this is the case, just skip to the next
        # underpopulated bin.
        if i in flatten(bins_to_merge):
            continue
        
        # Merge the under-populated bin with the successive bins until it has 'n_min' elements
        bins_to_merge_with_i = [i]    
        n = nn[i]
        j = i+1;
        while (n < n_min) and (j < n_bins):
            bins_to_merge_with_i.append(j)
            n += nn[j]
            j += 1
            
        bins_to_merge.append(bins_to_merge_with_i)
        bins_to_delete.append(bins_to_merge_with_i[1:])
        
        # Make vectors containing the properties of the bins to merge
        xx_i = xx.take(bins_to_merge_with_i)
        yy_i = yy.take(bins_to_merge_with_i)
        nn_i = nn.take(bins_to_merge_with_i)
        var_i = var.take(bins_to_merge_with_i);
            
        # Replace the value & position of the i-th bin with the average of the
        # To compute the value of the merged bin, we use a weighted average on
        # the number of pairs.
        xx[i] = xx_i.mean()                
        yy[i] = np.average( yy_i, weights=nn_i )
        
        # Compute the new variance.  This has two components:
        # 1) Contribution from the variances of the single bins:
        #        (nn^i)**2 var_i
        variance_contribution = np.dot(nn_i**2, var_i)
        # 2) contribution from the different averages of the different bins
        #    (it is the only one from a bin that contains only one data point):
        #        1/2 nn^i nn^j (<v>_i - <v>_j)**2
        delta_avg_squared = np.subtract.outer(yy_i,yy_i)**2
        average_contribution = 0.5 * np.einsum("i,j,ij", nn_i, nn_i, delta_avg_squared)
        # The variance also has a nice n**2 factor
        var[i] = 1./n**2 * ( variance_contribution + average_contribution )
        
        print "\n"
        print "Bins to merge                     " + str(np.array(bins_to_merge_with_i,dtype=np.int64))
        print "Number of data points             " + str(np.array(nn_i,dtype=np.int64))
        # print "Sigmas of bins                    " + str(np.sqrt(var_i))
        print "Sigma from the variances          " + str(np.sqrt(1./n**2 * variance_contribution))
        print "Sigma from the averages           " + str(np.sqrt(1./n**2 * average_contribution))
        print "New sigma                         " + str(np.sqrt(var[i]))
        
        # Update the number of data points in the merged bin. This line
        # should go last, after the line where we replace yy[i] and var[i]
        nn[i] = sum(nn.take(bins_to_merge_with_i))
        
    # Delete the bins that have been merged
    nn = np.delete(nn, flatten(bins_to_delete))
    xx = np.delete(xx, flatten(bins_to_delete))
    yy = np.delete(yy, flatten(bins_to_delete))
    var = np.delete(var, flatten(bins_to_delete))
    
    # Check that the bins all have more than 'n_min' data points
    if not all(nn[:-1] >= n_min):
        raise BaseException("Some bins are still under-populated.  Something went wrong.")
    
    # The last bin may still be underpopulated, because there
    # is no "next bin" to merge it with. Here we should take care of it
    # by merging it with the previous one.
    if nn[-1] < n_min:
        pass
        
    # Compute the standard deviation and the Gaussian error
    sigma = np.sqrt(var)
    err = np.sqrt(var/nn)
    
    # Sort out what to return, and return it
    return_list = [xx, yy, nn, var, sigma, err]
    
    if edges_requested:
        edges = np.delete(edges, flatten(bins_to_delete))        
        return_list.append(edges)
        
    return tuple(return_list)


# ===========================
# = Generate re-binned file =
# ===========================
filename = sys.argv[1]
npairs_treshold = int(sys.argv[2])
if len(sys.argv) > 3:
    out_filename = sys.argv[3]
else:
    out_filename = re.sub('.dat', '_rebin' + str(npairs_treshold) + '.dat', filename)

# Columns indices in the file
edges_col = 0
midpts_col = 1
npairs_col = 2
avgV_col = 3
varV_col = 4
sigmaV_col = 5
n_columns = 6

# Load file
data = np.loadtxt(filename, skiprows=1, usecols=range(n_columns))

# Exclude bins with zero pairs, in order to avoid getting NaNs
empty_bin_idxs = (data[:,npairs_col]==0).nonzero()
data = np.delete(data, empty_bin_idxs, axis=0)

# Re-bin the data
data_backup = data.copy()
(midpts, avgV, npairs, varV, sigmaV, errV, edges) = rebin(
    xx=data[:,midpts_col], yy=data[:,avgV_col], nn=data[:,npairs_col], n_min=npairs_treshold, var=data[:,varV_col], edges=data[:,edges_col])

# Compute error bars for the old binning
errV_old = data[:,sigmaV_col]/np.sqrt(data[:,npairs_col])

# Plot before/after
plt.close('all')

plt.subplot(2,1,1)
plt.errorbar(data[:,midpts_col], data[:,avgV_col], linestyle='None', yerr=errV_old, label=os.path.basename(filename))
# plt.plot(data[:,midpts_col], data[:,avgV_col], label=os.path.basename(filename))
plt.xscale('log')
# plt.xlim([np.min(data[:,edges_col]), np.min([np.max(data[:,edges_col]),40.])])
plt.ylim([np.min(data[:,avgV_col])-50,np.max(data[:,avgV_col])+50])
# plt.legend()

plt.subplot(2,1,2)
plt.errorbar(midpts, avgV, linestyle='None', yerr=errV, label=os.path.basename(out_filename))
plt.xscale('log')
# plt.xlim([0.005,40])
plt.ylim([np.min(data[:,avgV_col])-50,np.max(data[:,avgV_col])+50])
# plt.legend()

plt.show()


# Get the first line of the file containing the names of the columns
# and modify it to make it DataGraph-friendly
with open(filename, 'r') as fp:
    legend = fp.readline()
    legend = re.sub("_|\(|\)|\#", "", legend)

# Write to output file    
with open(out_filename, 'w') as fp:
    fp.write(legend)
    for i in range(len(midpts)):
        line = "%30f %30f %30d %30f %30f %30f\n" % (edges[i], midpts[i], int(npairs[i]), avgV[i], varV[i], sigmaV[i])
        fp.write(line)



