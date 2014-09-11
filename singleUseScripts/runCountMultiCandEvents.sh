#!/bin/bash

root -b -x -q countMultiCandEvents.C+\(\"/raid/raid7/jhugon/hmmNtuples/stage2/ggHmumu125_8TeV.root\"\)
root -b -x -q countMultiCandEvents.C+\(\"/raid/raid7/jhugon/hmmNtuples/stage2/vbfHmumu125_8TeV.root\"\)
root -b -x -q countMultiCandEvents.C+\(\"/raid/raid7/jhugon/hmmNtuples/stage2/wHmumu125_8TeV.root\"\)
root -b -x -q countMultiCandEvents.C+\(\"/raid/raid7/jhugon/hmmNtuples/stage2/zHmumu125_8TeV.root\"\)

root -b -x -q countMultiCandEvents.C+\(\"/cms/data/store/user/jhugon/hmumu/stage2/DYJetsToLL_8TeV.root\"\)
root -b -x -q countMultiCandEvents.C+\(\"/cms/data/store/user/jhugon/hmumu/stage2/ttbar_8TeV.root\"\)
root -b -x -q countMultiCandEvents.C+\(\"/cms/data/store/user/jhugon/hmumu/stage2/WW_8TeV.root\"\)
root -b -x -q countMultiCandEvents.C+\(\"/cms/data/store/user/jhugon/hmumu/stage2/WZ_8TeV.root\"\)
root -b -x -q countMultiCandEvents.C+\(\"/cms/data/store/user/jhugon/hmumu/stage2/ZZ_8TeV.root\"\)

root -b -x -q countMultiCandEvents.C+\(\"/cms/data/store/user/jhugon/hmumu/stage2/SingleMuRun2012Av1-22Jan2013.root\"\)
root -b -x -q countMultiCandEvents.C+\(\"/cms/data/store/user/jhugon/hmumu/stage2/SingleMuRun2012Bv1-22Jan2013.root\"\)
root -b -x -q countMultiCandEvents.C+\(\"/cms/data/store/user/jhugon/hmumu/stage2/SingleMuRun2012Cv1-22Jan2013.root\"\)
root -b -x -q countMultiCandEvents.C+\(\"/cms/data/store/user/jhugon/hmumu/stage2/SingleMuRun2012Dv1-22Jan2013.root\"\)
