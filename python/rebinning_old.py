#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
import re

def flatten( list ):
    """Turn a list of lists into a list with the elements of the sublists.
    Trick taken from http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
    """
    return sum(list,[])


def rebin( xx, yy, nn, n_min, edges=None, var=None, sigma=None, err=None ):
    """Rebin an array imposing a minimum number of points in each bin.
    
    Input:
    xx, yy, nn -> 1D array containing, respectively, the bin-midpoints, the values of the 
        bins and the number of data points for each bin.
           
    n_min -> minimum number of data points each bin must contain in the returned arrays.
    
    edges -> 1D array with the edges associated to the original binning.  If present, the
        function will return also the re-binned edges.

    var -> 1D array with the variance associated to each bin.  If present, the function will return
        also the re-binned variance.

    sigma -> 1D array with the standard deviation associated to each bin.  If present, the function will return
        also the re-binned standard deviation.

    err -> 1D array with the uncertainty associated to each bin.  If present, the function will return
        also the re-binned uncertainty, assuming statistical errors only.
    
    Output:
    (xx, yy, nn, [edges, variance, sigma, err])
    
    TODO: When we merge two bins with, respectively, a population of 'n1' and 'n2' data points,
        we compute the variance of the resulting bin (of N=n1+n2 points) using the formula:

            var = 1/(n1+n2) (var1*n1 + var2*n2)
        
        which leads to an error estimate of:
        
            err = sqrt(var/N) = 1/(n1+n2) sqrt( (err1*n1)^2 + (err2*n2)^2 ) .
            
        This seems ok, but when one of the original bins has just one point, there is a paradox.
        What happens is that if n2=1, then <z^2> = <z>^2 = z^2 and therefore var2=0.
        As a result, we have that the total variance is 1/(n1+1)*var1*n1, which means that
        the single point in the second bin does not contribute to the variance of the merged bin.
        
        We need to fix this.  An idea is to compute the variance with the original formula for the
        merged bin.  Let's call x, y the two points in the first bin, and z the only point in the
        second bin.
        
            var = 1/(n1+1) * (x^2 + y^2 + z^2) - 1/(n1+1)^2 (x + y + z)^2

        If we substitute x^2 + y^2 = n1*var1,  x + y = n1*avg1,  and z = avg2,  we have an expression
        for the the variance of the resulting bin in terms of var, avg1, n1, avg2 which are known
        quantities.

    """

	# Check that we have enough data points so that can make at least one bin
    assert np.sum(nn) > n_min
    n_total = np.sum(nn)    

	# Check that the input is consistent
    assert len(xx) == len(yy) == len(nn)
    n_bins = len(xx)
    
    # Understand whether the user wants to compute the edges, variances, and errors on 
    # the average, with respect to the new bins.
    edges_requested = isinstance(edges, list) or isinstance(edges, np.ndarray)
    var_requested = isinstance(var, list) or isinstance(var, np.ndarray)
    sigma_requested = isinstance(sigma, list) or isinstance(sigma, np.ndarray)    
    err_requested = isinstance(err, list) or isinstance(err, np.ndarray)    
    
    # Work with arrays
    xx = np.asarray(xx, dtype=np.float64)
    yy = np.asarray(yy, dtype=np.float64)
    nn = np.asarray(nn, dtype=np.float64)
    if edges_requested:
        edges = np.asarray(edges, dtype=np.float64)
    if var_requested:
        var = np.asarray(var, dtype=np.float64)
    if sigma_requested:
        sigma = np.asarray(sigma, dtype=np.float64)
    if err_requested:
        err = np.asarray(err, dtype=np.float64)

    
    
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
            
        # Replace the i-th bin with the merged bin, and remove all the bins from i+1 to j.
        # To compute the 'y' value in the bin, we use a weighted average on the number of pairs.
        xx[i] = xx.take(bins_to_merge_with_i).mean()                
        yy[i] = np.average( yy.take(bins_to_merge_with_i), weights=nn.take(bins_to_merge_with_i) )
        if var_requested:
            var[i] = np.average( var.take(bins_to_merge_with_i), weights=nn.take(bins_to_merge_with_i) )
        if sigma_requested:
            sigma[i] = np.sqrt(np.dot ( sigma.take(bins_to_merge_with_i)**2, nn.take(bins_to_merge_with_i) ) / sum(nn.take(bins_to_merge_with_i)) )    
        if err_requested:
            err[i] = np.sqrt(np.dot ( err.take(bins_to_merge_with_i)**2, nn.take(bins_to_merge_with_i)**2 )) / sum(nn.take(bins_to_merge_with_i))  
            

        print "\n"
        print "Bins to merge                     " + str(np.array(bins_to_merge_with_i,dtype=np.int64))
        print "Number of data points             " + str(np.array(nn.take(bins_to_merge_with_i),dtype=np.int64))
        # print "Sigmas of bins                    " + str(np.sqrt(var_i))
        print "New sigma                         " + str(np.sqrt(var[i]))

            
        # This line should go last, after the line where we replace yy[i]
        nn[i] = sum(nn.take(bins_to_merge_with_i))

    # Delete the bins that have been merged
    xx = np.delete(xx, flatten(bins_to_delete))
    yy = np.delete(yy, flatten(bins_to_delete))
    nn = np.delete(nn, flatten(bins_to_delete))
    
    # Check that the bins all have more than 'n_min' data points
    if not all(nn[:-1] >= n_min):
        raise BaseException("Some bins are still under-populated.  Something went wrong.")
    
    # The last bin may still be underpopulated, because there
    # is no "next bin" to merge it with. Here we should take care of it
    # by merging it with the previous one.
    if nn[-1] < n_min:
        pass

    # Sort out what to return, and return it
    return_list = [xx,yy,nn]

    if edges_requested:
        edges = np.delete(edges, flatten(bins_to_delete))        
        return_list.append(edges)
    if var_requested:
        var = np.delete(var, flatten(bins_to_delete))                
        return_list.append(var)
    if sigma_requested:
        sigma = np.delete(sigma, flatten(bins_to_delete))                
        return_list.append(sigma)
    if err_requested:
        err = np.delete(err, flatten(bins_to_delete))                
        return_list.append(err)
    
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
(midpts, avgV, npairs, edges, varV, sigmaV) = rebin(
    xx=data[:,midpts_col], yy=data[:,avgV_col], nn=data[:,npairs_col], n_min=npairs_treshold,
    edges=data[:,edges_col], var=data[:,varV_col], sigma=data[:,sigmaV_col])

# Compute error bars
errV_old = data[:,sigmaV_col]/np.sqrt(data[:,npairs_col])
errV = sigmaV/np.sqrt(npairs)

# Plot before/after
plt.subplot('211')
plt.errorbar(data[:,midpts_col], data[:,avgV_col], linestyle='None', yerr=errV_old)
plt.xscale('log')
plt.xlim([0.005,40])
plt.ylim([np.min(data[:,avgV_col])-50,np.max(data[:,avgV_col])+50])

plt.subplot('212')
plt.errorbar(midpts, avgV, linestyle='None', yerr=errV)
plt.xscale('log')
plt.xlim([0.005,40])
plt.ylim([np.min(data[:,avgV_col])-50,np.max(data[:,avgV_col])+50])

plt.show()


# Get the first line of the file containing the names of the columns
with open(filename, 'r') as fp:
    legend = fp.readline()

# Write to output file    
with open(out_filename, 'w') as fp:
    fp.write( re.sub('\#', ' ', legend) )
    for i in range(len(midpts)):
        line = "%30f %30f %30d %30f %30f %30f\n" % (edges[i], midpts[i], int(npairs[i]), avgV[i], varV[i], sigmaV[i])
        fp.write(line)


# x=[-1,0,1,2,3,4,5,6];
# edges=[-1.5,-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5];
# y=[-20,-15,-10,-5,0,5,10,15];
# n=[30,30,10,20,30,100,15,30];
# err=[1,2,3,4,5,6,7,8]
# b=rebin(x,y,n,50,edges=edges,err=err)
# print b


