#!/bin/bash
#
# get_millennium.sh
# Created by Alicia Bueno Belloso and Guido W. Pettinari
# Last modified by Guido W. Pettinari on 22/01/2012
# 
# Download from a remote server galaxy catalogues from the Millennium simulation.
# usage:  get_millennium <model> <output_folder> [<start_from_id>]
#
# The main output file will be stored in 'output_folder' and will be named according
# to the chosen 'model' and to the cuts defined inside this file in the 'cuts'
# variable.  The criteria of the performed SQL query (cuts, columns, table)
# are hard-coded in this script just after the parsing of the arguments.
#  
# The main feature of this script is that it automatically resumes any timed-out
# query.  This is particularly useful when downloading the guo2010a catalog from
# the Garching mirror, since it only allows 7-minute queries.
# The username and password to connect to the Millennium remote servers
# are hard coded in the function below, query_millennium.
#
# The following options for the 'model' argument are accepted:
# font2008a (see http://arxiv.org/abs/0807.0001)
# guo2010a (see http://arxiv.org/abs/1006.0106)
#
# If the 'start_from_id' argument is provided, then the query will only consider
# galaxies where galaxy_id > $start_from_id.  This is useful to resume a previous
# run of the script that was interrupted for any reason.


# Support function for the script.
# Submit a query to a custom server of the Millennium simulation.
# The supported servers are the Garching one (garching) and the
# Durham mirror (durham)
function query_millennium () {

  query="$1"
  server="$2"

  if [[ $server == 'durham' ]]; then
    ### Define arguments for wget (Durham)
    url='http://galaxy-catalogue.dur.ac.uk:8080/MyMillennium?action=doQuery&SQL='
    action='wget'
    ## My login
    user='--http-user=pettinari'
    password='--http-passwd=4btH9r62'
    ## Niko's login
    # user='--http-user=meures'
    # password='--http-passwd=MGMduz05'

  elif [[ $server == 'garching' ]]; then
    
    ### Define arguments for wget (Garching)
    url='http://gavo.mpa-garching.mpg.de/MyMillennium3?action=doQuery&SQL='
    action='wget'
    user='--http-user=abueno'
    password='--http-passwd=tJPrUjuVE'
  
  else
    echo "The 'server' argument must be either 'durham' or 'garching', not $server." >&2
    exit 1
  fi
  
  # When building the query, make sure to correctly encode white-spaces as %20
  query=${query// /%20}

  # Build and execute command
  args="-O - $user $password $url$query"
  cmd="$action $args"
  echo "*** $cmd" >&2
  $cmd
  
}




### Parse arguments
if [[ $# < 2 ]]; then
  echo "usage:  $0 <model> <output_folder> [<start_from_id>]"
  exit 1
fi
model="$1"
output_folder="$2"
start_from_id="$3"

if ! [[ -d $output_folder ]]; then
  printf "Directory doesn't exist, creating directory %s\n" "$output_folder"
  mkdir "$output_folder"
fi





# ===================================
# = Modify here to change the query =
# ===================================

if [[ $model == 'font2008a' ]]; then
  ### Define arguments for the SQL query (Font2008)
  server='durham'
  table='dgalaxies..font2008a'
  columns="galaxyid,x,y,z,vx,vy,vz,type,mhalo,stellarMass,r_SDSS"
  # mass_cut='mhalo>73e10'
elif [[ $model == 'guo2010a' ]]; then
  ### Define arguments for the SQL query (Guo2010a)
  server='garching'
  table='guo2010a..mr'
  columns="galaxyid,x,y,z,velX,velY,velZ,type,np,mvir,stellarMass,r_mag,rDust"
  # columns="galaxyid,type,np,mvir,haloid,r_mag"
  # columns="galaxyid,x,y,z,velx,vely,velz,centralmvir,mvir,r_mag"
  # mass_cut='centralmvir>73'
else
  echo "The 'model' argument must be either 'font2008a' or 'guo2010a', not $model."
  exit 1
fi

### Cuts

# For redshift z=0, use snapnum=63. For redshift z=0.989, use snapnum=41.
# For redshift z=1.504, use snapnum=36. Fore z=0.5085 use snapnum=48.
# For a table of equivalences
# between redshift and snapshots (see select snapnum,redshift from MField..Snapshots)
snapnum='41'

# Here is the place to add custom cuts or modify the existing ones.
# Memorize in 'cuts' all the cuts that will be made in the the SQL query.
# The results file will be named using the content of 'cuts'.
cuts="snapnum=$snapnum"
# cuts="$mass_cut and snapnum=$snapnum"




# =============================
# = Do not modify from now on =
# =============================

### Query building

# Function that generates the query.  When given an argument, the query will
# contain a cut on the galaxyid in order to match the last query
function generate_query () {
  if [[ $1 == '' ]]; then
    echo "select $columns from $table where $cuts"
  else
    echo "select $columns from $table where $cuts and galaxyid>$1"
  fi
}

# The first query won't have any cut in the galaxyid unless it is specified
# from the command line, in which case $start_from_id would be different from ''
query=`generate_query $start_from_id`
printf "The first query to be run will be:\n     ->   %s\n" "$query"

# Function to generate the file names where the intermediate pieces of
# the simulation will be stored
function temp_root() {
  echo "$output_folder/temp$1.dat"
}


# ==============
# = Main cycle =
# ==============
# When finish==1, the main cycle will stop
finish=0
counter=1

while [[ $finish -eq 0 ]]; do

  temp_file=`temp_root $counter`
  printf "*** w-getting file %s...\n" "$temp_file"
  # printf "The new query is %s\n" "$query"  
  query_millennium "$query" "$server" > "$temp_file"
  printf "Done!\n"
  
  # Check if the query did time out.  If yes, rerun it starting from the right
  # point, otherwise we are done
  last_line=`tail -n1 $temp_file`
  printf "*** Last line of file is %s\n" "$last_line"
  echo;
  
  if [[ $last_line == *TIMED\ OUT* || $last_line == *SQLEXCEPTION* ]]; then

    printf "*** Previous query timed out.  Be calm, and carry on with another query.\n"
    # Use awk to get the first field of the last line
    last_galaxy_id=`sed '/^#/d' "$temp_file" | awk -F"," 'END{print $1}'`
    printf "*** Last galaxy id is %s\n" "$last_galaxy_id"
    # Update the query, so that the next one starts from the last galaxyid we got
    query=`generate_query $last_galaxy_id`
    printf "    The new query is %s\n" "$query"  
    # Update the counter
    counter=$(($counter+1))
   
  else 
    finish=1
  fi
  
done

# Debug
printf "*** Out of the main loop.\n"
last_line=`tail -n1 $temp_file`
printf "*** Last line of last file is %s\n" "$last_line"
printf "*** Line count of output directory BEFORE processing:\n"
wc -l $output_folder/*
echo

# ===================================
# = Combine files into one big file =
# ===================================
# Generate a unique file from the intermediate ones created in the main loop.
# This unique file contains only data, no comment lines.
# IMPORTANT: If you want to run only this part of the code, just comment the while loop above,
# and set 'counter' equal to the number of temp files you've got.
echo "**** START PROCESSING FILES"

# Build filename of final file
output_path="$output_folder/galaxies_${model}_${columns//","/_}_where_${cuts//" and "/_}.dat"
output_path=${output_path//"="/}
# Since we shall append the results to the file, we issue a warning ig
# it already exists
printf "*** Final result will be stored in %s\n" "$output_path"
if [[ -f $output_path ]]; then
  mv "$output_path" "$output_path.backup"
  printf "    WARNING: file already exists.  The existing file was renamed as %s.\n" "$output_path.backup"
fi

for (( i = 1; i <= $counter; i++ )); do

  temp_file=`temp_root $i`
  # Append the cleaned file to the final file, removing non-data lines
  printf "    Appending file %s... " "$temp_file"
  cat "$temp_file" | sed '/^[^0-9]/d' >> "$output_path"
  # If speed is an issue and you're sure that the non-data lines all
  # start with a hash, then use the following line instead:
  # cat "$temp_file" | sed '/^#/d' >> "$output_path"  
  echo "done."
  
done

# Debug
printf "*** The final file '%s' is now complete.\n" "$output_path"
printf "*** Line count of output directory AFTER processing:\n"
wc -l $output_folder/*






