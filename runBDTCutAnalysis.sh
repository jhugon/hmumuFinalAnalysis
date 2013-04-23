#!/bin/bash

rm -f shapes/*
rm -f statsCards/*
rm -f statsInput/*
rm -f statsOutput/*

nice ./makeCards.py --cutOpt
cd statsCards/
nice bash notlxbatch.sh
cd ..
cp statsCards/*.out statsInput/.
nice ./makeLimitPlots.py --cutOpt
