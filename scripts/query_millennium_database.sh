#!/bin/bash
#
# Submit a query to a custom server of the Millennium simulation.
# usage:    query_millennium_database <query> <server>
#
# The supported servers are the Garching main server (garching) and the
# Durham mirror (durham)

### Parse arguments
if [ $# != 2 ]; then
  echo "usage:  $0 <query> <server>"
  exit 1
fi
query="$1"
server="$2"

if [[ $server == 'durham' ]]; then
  ### Define arguments for wget (Durham)
  url='http://galaxy-catalogue.dur.ac.uk:8080/MyMillennium?action=doQuery&SQL='
  action='wget'
  user='--http-user=meures'
  password='--http-passwd=MGMduz05'

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
echo $cmd >&2
$cmd