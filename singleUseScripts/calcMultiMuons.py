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

def readLogFile(sampleName):
  result = None
  if sampleName == "data":
    runs = [
                "SingleMuRun2012Av1",
                "SingleMuRun2012Bv1",
                "SingleMuRun2012Cv1",
                "SingleMuRun2012Dv1",
            ]
    for run in runs:
        runData = readLogFile(run)
        if not result:
          result = runData
        else:
          for cat in result.keys():
            for i in [0,1]:
              result[cat][i] += runData[cat][i]
  else:
    result = {}
    if sampleName == "TT":
      sampleName = "ttbar"
    files = list(glob.glob(sampleName+"*.txt"))
    assert(len(files)==1)
    f = open(files[0])
    for line in f:
      #print line,
      match = re.match(r"^(\w+)\:\s+(\d+)\s+(\d+)$",line)
      assert(match)
      result[match.group(1)] = [int(match.group(2)),int(match.group(3))]
    f.close()
  return result
    

if __name__ == "__main__":
  categories = [
                    "MuonSelected",
                    "01JetTightBB",
                    "01JetTightBO",
                    "2JetVBFTight",
                    "2JetGFTight",
               ]
  samples = xsecDict.keys()+["data"]
  samples.sort()
  data = {}
  for sample in samples:
    data[sample] = readLogFile(sample)

  samples.remove("data")
  samples2 = ["SigTotal","BkgTotal","data"]
  data["SigTotal"] = {}
  data["BkgTotal"] = {}
  for cat in categories:
    data["SigTotal"][cat] = [0.,0.]
    data["BkgTotal"][cat] = [0.,0.]

  for i in [0,1]:
    for sample in samples:
      sf = luminosity*1000.*xsecDict[sample]/originalNEventsDict[sample]
      for cat in categories:
        if len(sample)==2:
          data["BkgTotal"][cat][i] += data[sample][cat][i]*sf
        else:
          data["SigTotal"][cat][i] += data[sample][cat][i]*sf

  for i,trial in zip([0,1],["Number of Independent Events","Number of >1 Higgs Candidate Events"]):
    print "="*6*15
    print trial+" for 19.7 fb^-1 at 8 TeV"
    print "="*6*15
    outStr = "{0:15}".format("")
    for cat in categories:
      outStr += "{0:>15}".format(cat)
    print outStr
    for sample in samples:
      outStr = "{0:15}".format(sample)
      sf = luminosity*1000.*xsecDict[sample]/originalNEventsDict[sample]
      for cat in categories:
        outStr += "{0:15.3g}".format(data[sample][cat][i]*sf)
      print outStr
    print "-"*6*15
    for sample in samples2:
      outStr = "{0:15}".format(sample)
      for cat in categories:
        outStr += "{0:15.3g}".format(data[sample][cat][i])
      print outStr
    print "="*6*15
    print
