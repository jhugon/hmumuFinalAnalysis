#!/bin/bash

hn=$(hostname)
un=$(whoami)
if [[ "$hn" == *cern* ]]; then
  njobs=`bjobs | wc -l`
  njobsRunning=`bjobs | grep "RUN" | wc -l`
  ((njobs -= 1))
else
  rawJobOut=$(qstat -u $un -t)
  njobsQ=$(echo "$rawJobOut" | grep '^[0-9].* Q ' | wc -l)
  njobsRunning=$(echo "$rawJobOut" | grep '^[0-9].* R ' | wc -l)
  njobs=""
  ((njobs = njobsQ + njobsRunning))
fi
echo "N jobs:             "$njobs
echo "N jobs Running:     "$njobsRunning
echo "Running Percentage: "$(python -c "print float($njobsRunning)/$njobs")
