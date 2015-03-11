#!/bin/bash

## Connect to a remote ssh server and run ./pair_selector.
## You need to have defined a script that logs into the ssh server
## by just executing the name of the server.  This script is supposed
## to be located in /Users/coccoinomane/my_projects/shell_scripts
##
## IMPORTANT: make sure that the remote params.ini file was updated on
## with the latest version!

if [ $# != 3 ]; then
  echo "usage:  $0 <server> <params.ini file> <n_threads>"
  exit 1
fi
server=$1
params_file=$2
n_threads=$3

## Remote home & project directories
RHOME=/Users/pettinag
remote_project_dir="$RHOME/my_projects/alcock_paczynski_working_copy/Code"
pair_selector_cmd="./pair_selector $params_file"

## Local path of the SSH script; must execute its argument into remote machine,
## possibly without asking for password
ssh_script=/Users/coccoinomane/my_projects/shell_scripts/$server
log_file="${server}_remote_script.log"

## Commands to be executed remotely
ssh_cmd="\
  cd $remote_project_dir;\
  export OMP_NUM_THREADS=$n_threads;\
  time $pair_selector_cmd >& $log_file\
"
## We are going to execute the following list of commands:
echo "We are going to execute remotely the following commands:"
echo $ssh_cmd

## Execute the script
$ssh_script "$ssh_cmd"