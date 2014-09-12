#!/usr/bin/env python

import glob
import ROOT as root
import singleHelpers
from helpers import *

root.gROOT.SetBatch(True)


#####################################
## Configuration (all for 8TeV)

# Scale factor for higgs cross-sections
SIGSF=1.

luminosity = 19.712 # fb^-1 Pixel-Lumi for 22Jan2013 ReRECO 

originalNEventsDict = {
  # Signal Needs to be checked!!
  "ggHmumu":100000,
  "vbfHmumu":100000,
  "wHmumu":100000,
  "zHmumu":100000,

  "DY":30086987,
  "TT":6921652,

  "WW":5218045,
  "WZ":10000283,
  "ZZ":10320000,
}

# in pb
xsecDict = {
  "ggHmumu":0.0042944*SIGSF,
  "vbfHmumu":0.00034716*SIGSF,
  "wHmumu":0.000153252*SIGSF,
  "zHmumu":8.6746e-05*SIGSF,

  "DY":3503.71,
  "TT":225.197,

  "WW":54.838,
  "WZ":33.21,
  "ZZ":17.654,
}

multiMuonCountsDict = {
  "ggHmumu":0,
  "vbfHmumu":0,
  "wHmumu":1929,
  "zHmumu":1143,

  "DY":0,
  "TT":2,

  "WW":0,
  "WZ":478,
  "ZZ":625,
}

muonSelCountsDict = {
  "ggHmumu":50987,
  "vbfHmumu":55610,
  "wHmumu":47940,
  "zHmumu":48849,

  "DY":46567,
  "TT":9100,

  "WW":4989,
  "WZ":6418,
  "ZZ":9247,
}

if __name__ == "__main__":
  totalEvents = 0.
  for sampleName in multiMuonCountsDict:
    toprint = "{0:10}".format(sampleName)
    #count = multiMuonCountsDict[sampleName]
    count = muonSelCountsDict[sampleName]
    toprint += "{0:10}".format(count)
    count *= luminosity*1000.*xsecDict[sampleName]/originalNEventsDict[sampleName]
    toprint += "{0:20f}".format(count)
    print toprint
    totalEvents += count
  print totalEvents
