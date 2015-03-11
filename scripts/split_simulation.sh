### "usage:  split_file.sh <grid_size> <n_bins>"
###
### Bash script to split the Millenium simulation into 1000 sub-files using
### the 'file_splitter' utility.
###
### This script will create n_bins^3 files. If n_bins > 6, most likey you will
### need to increase the maximum number of open-able files with the command
### ulimit -S -n n_files
### where n_files should be at least n_bins^3.


## Parse arguments
if [ $# != 2 ]; then
  echo "usage:  $0 <grid_size> <n_bins>"
  echo "N.B.  If you get segmentation fault, most likely you need to"
  echo "      increment the maximum number of openable files with"
  echo "      ulimit -S -n n_files"
  echo "      where n_files should be at least n_bins^3."
  exit 1
fi
grid_size=$1
n_bins=$2

## Here we update the limit imposed 
## by the shell on the maximum number of files to n_bins^3.
## NOT WORKING!
# n_files=$((n_bins*n_bins*n_bins))
# if [ $n_files -gt 127 ]; then
#   ulimit -S -n $n_files
# fi

## Stuff you may want to change
in_dir="../Data"
out_dir="../Data/data_blocks_z2"
in_file="$in_dir/galaxies_font2008_z_1504.dat"
working_dir="."
bin="./file_splitter"
out_files_format="$out_dir/gal_blocks_%d_%d_%d.dat"

## Build the command
cmd="$bin $in_file $out_files_format $grid_size $n_bins"
echo $cmd

## Execute the command
$cmd