#!/bin/bash
#
# Strip the first column from a file.  Useful to remove the galaxyid column
# from the files obtained through "get_millennium.sh"

if [[ $# != 2 ]]; then
  echo "usage:  $0 <input_file> <output_file>"
  exit 1
fi
input_file=$1
output_file=$2

cut -d, -f2- $input_file > $output_file