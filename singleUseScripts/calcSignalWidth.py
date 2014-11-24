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

def getMinWidth(x,prob=0.683):
  """
  Find the width of the minimum interval containing "prob" fraction
  of values of the array "x"
  """
  xTmp = sorted(x)
  xLen = len(xTmp)
  widthList = []
  nPoints = int(numpy.floor(prob*xLen))
  #print "Points: ",nPoints
  for i in range(xLen-nPoints+1):
    width = xTmp[i+nPoints-1]-xTmp[i]
    #print i,xTmp[i+nPoints],xTmp[i],width
    #print i,i+nPoints-1,xTmp[i],xTmp[i+nPoints-1],width
    widthList.append(width)
  return numpy.amin(widthList)

inputFileNameGF = getDataStage2Directory()+"/ggHmumu125_8TeV.root"
treeGF = root.TChain("outtree")
treeGF.AddFile(inputFileNameGF)
treeGF.SetCacheSize(10000000);
treeGF.AddBranchToCache("*");

inputFileNameVBF = getDataStage2Directory()+"/vbfHmumu125_8TeV.root"
treeVBF = root.TChain("outtree")
treeVBF.AddFile(inputFileNameVBF)
treeVBF.SetCacheSize(10000000);
treeVBF.AddBranchToCache("*");

analyses = []
jet2PtCuts = " && jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40."
jet01PtCuts = " && !(jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40.)"
categoriesAll = ["BB","BO","BE","OO","OE","EE"]
analyses += [["Jets01PassPtG10"+x,  "dimuonPt>10." +jet01PtCuts] for x in categoriesAll]
analyses += [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]
analyses += [["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts]]
analyses += [["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
analyses += [["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]

print "{0:30} {1:4} {2:4}".format("","GF","VBF")
for analysis in analyses:
  analysisName = analysis[0]
  selectionString = analysis[1]
  fullSelectionString = treeCut(analysisName,selectionString,eventWeights=False,muonRequirements=True)
  massListGF = getListOfVarFromTree(  treeGF,
                                    "dimuonMass",
                                    fullSelectionString
                                 )
  massListVBF = getListOfVarFromTree(  treeVBF,
                                    "dimuonMass",
                                    fullSelectionString
                                 )
  relWidthGF = 0.5*getMinWidth(massListGF)/125.0
  relWidthVBF = 0.5*getMinWidth(massListVBF)/125.0
  print "{0:30} {1:4.1%} {2:4.1%}".format(analysisName,relWidthGF,relWidthVBF)













