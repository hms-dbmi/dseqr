#!/bin/bash

cids=`cat cids.csv`
cidstot=`echo "$cids" | wc -w`

echo "have $cidstot cids"

cidnum=0
t0=`date +%s%N | cut -b1-13`

for cid in $cids
do
  # download if doesn't exist
  outf="views/$cid.json"
  if [ ! -f $outf ]
  then
    curl -s https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/$cid/JSON/ > views/$cid.json
  fi


  ((cidnum++))
  echo "Finished $cid ($cidnum of $cidstot)"

  # make sure only downloading 5 per second
  if [ $(( cidnum % 5 )) == 0 ]
  then
    # wait remaining up to 1 second
    t1=`date +%s%N | cut -b1-13`
    elapsed=$(($t1-$t0))
    remains=$((1000-$elapsed))
    if [ $remains -gt 0 ]
    then
      sleep_seconds=`bc -l <<< $remains/1000`
      echo sleeping $sleep_seconds
      sleep $sleep_seconds
    fi

    # reset start time
    t0=`date +%s%N | cut -b1-13`
  fi
done
