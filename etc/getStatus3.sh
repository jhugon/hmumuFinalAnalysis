#!/bin/bash

njobs=`bjobs | wc -l`
njobsRunning=`bjobs | grep "RUN" | wc -l`
((njobs -= 1))
echo "N jobs:             "$njobs
echo "N jobs Running:     "$njobsRunning
echo "Running Percentage: "$(python -c "print float($njobsRunning)/$njobs")
