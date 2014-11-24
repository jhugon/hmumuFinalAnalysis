#!/usr/bin/env python

import singleHelpers
from helpers import *
import makeCards
import random
import sys

def getListOfVarFromTree(tree,varName,cutString=""):
#  print cutString
  result = []
  entryListName = "myEntryList{0:08x}".format(random.getrandbits(32))
  #print entryListName
  tree.Draw(">>"+entryListName,cutString,"entryList")
  entryList = root.gDirectory.Get(entryListName)
  for iEvent in range(tree.GetEntries()):
    if entryList.Contains(iEvent):
      tree.GetEntry(iEvent)
      result.append(getattr(tree,varName))
  return result

def getQuantileWidth(data,sigmas=1.):
  quantiles = numpy.percentile(massList,
                      [
                          100.*scipy.stats.norm.cdf(-sigmas),
                          50.,
                          100.*scipy.stats.norm.cdf(sigmas)
                      ]
                    )
  #print quantiles, quantiles[2]-quantiles[0]
  return (quantiles[2]-quantiles[0])/2.

inputFileName = getDataStage2Directory()+"/ggHmumu125_8TeV.root"
tree = root.TChain("outtree")
tree.AddFile(inputFileName)
tree.SetCacheSize(10000000);
tree.AddBranchToCache("*");

analyses = []
jet2PtCuts = " && jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40."
jet01PtCuts = " && !(jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40.)"
categoriesAll = ["BB","BO","BE","OO","OE","EE"]
analyses += [["Jets01PassPtG10"+x,  "dimuonPt>10." +jet01PtCuts] for x in categoriesAll]
analyses += [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]
analyses += [["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts]]
analyses += [["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
analyses += [["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]

for analysis in analyses:
  analysisName = analysis[0]
  selectionString = analysis[1]
  fullSelectionString = treeCut(analysisName,selectionString,eventWeights=False,muonRequirements=True)
  massList = getListOfVarFromTree(  tree,
                                    "dimuonMass",
                                    fullSelectionString
                                 )
  print "{0:30} {1:4.1%}".format(analysisName,getQuantileWidth(massList)/125.0)













