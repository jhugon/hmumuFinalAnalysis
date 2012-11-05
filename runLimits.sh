#!/bin/bash

rm -f shapes/*
rm -f statsCards/*
rm -f statsInput/*
rm -f statsOutput/*

./makeCards.py
cd statsCards/
bash notlxbatch.sh
cd ..
cp statsCards/*.out statsInput/.
./makeShapePlots.py
./makeLimitPlots.py
