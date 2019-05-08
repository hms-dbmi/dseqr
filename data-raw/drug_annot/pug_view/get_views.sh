#!/bin/bash

cids=`cat cids.csv`
cidstot=`echo "$cids" | wc -w`

echo "have $cidstot cids"

cidnum=0

for cid in $cids
do
  # download if doesn't exist
  outf="views/$cid.json"
  if [ ! -f $outf ]
  then
    curl -s https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/$cid/JSON/ > views/$cid.json
  else
    ((cidnum++))
    echo "Finished $cid ($cidnum of $cidstot)"
    continue
  fi

  ((cidnum++))
  echo "Finished $cid ($cidnum of $cidstot)"

  # wait 0.2 seconds so that only 5 per second
  sleep 0.2
done
