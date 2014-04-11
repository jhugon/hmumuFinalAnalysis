#!/usr/bin/env python

import singleHelpers
from helpers import *
from xsec import *
import cPickle
import sys
import os
import os.path
import glob
import re
import numpy
from numpy import mean,median

from fitBiasStudy import mergeDicts

from biasPklToMu import getSMSigCounts

errorsDictFile = open("pklfiles/oneSig.pkl")
ERRORSDICT = cPickle.load(errorsDictFile)
errorsDictFile.close()

def getBiasStudies(globStr):
  """
  loads in biasstudy raw .pkl files
  returns a list of combined dictionaries
  """
  result = []
  inputPklFiles = glob.glob(globStr)
  tmpInputPklFiles = inputPklFiles
  foundJobGroupPkl = False
  foundNotJobGroupPkl = False
  for inputPkl in inputPklFiles:
    if "_jobGrp" in inputPkl:
      foundJobGroupPkl = True
    else:
      foundNotJobGroupPkl = True
  if foundNotJobGroupPkl and foundJobGroupPkl:
      print "Error: found .pkl files containing '_jobGrp' and not containing '_jobGrp'."
      print "  Can only process one or the other.  "
      print "  Please delete one or the other and try again. Exiting."
      sys.exit(1)
  if foundNotJobGroupPkl:
    for inputPkl in inputPklFiles:
      print "Running over input pkl file: "+inputPkl
      result.append(inputPkl)
  else:
    # Identify basenames to combine job groups
    basenames = set()
    for inputPklFn in inputPklFiles:
      match = re.match(r"(.+jobGrp)[\d]+\.pkl",inputPklFn)
      assert(match)
      tmpBase = match.group(1)
      if not (tmpBase+"*.pkl") in basenames:
          basenames.add((tmpBase+"*.pkl"))
    for globStr in basenames:
      fns = glob.glob(globStr)
      resultData = None
      for tmpFn in fns:
        tmpF = open(tmpFn)
        tmpD = cPickle.load(tmpF)
        if resultData == None:
          resultData = tmpD 
        else:
          mergeDicts(resultData,tmpD,True)
        tmpF.close()
      result.append(resultData)
  # Sort result by category name
  orderDef = [
    "CombSplitAll",
    "Jets01SplitCatAll",
    "Jet2SplitCutsGFSplit",
    "Jets01PassCatAll" ,
    "Jets01FailCatAll" ,

    "Jets01PassPtG10BB",
    "Jets01PassPtG10BO",
    "Jets01PassPtG10BE",
    "Jets01PassPtG10OO",
    "Jets01PassPtG10OE",
    "Jets01PassPtG10EE",
                          
    "Jets01FailPtG10BB",
    "Jets01FailPtG10BO",
    "Jets01FailPtG10BE",
    "Jets01FailPtG10OO",
    "Jets01FailPtG10OE",
    "Jets01FailPtG10EE",

    "Jet2CutsVBFPass",
    "Jet2CutsGFPass",
    "Jet2CutsFailVBFGF",
  ]
  result.sort(key=lambda x: orderDef.index(x['meta']['catName']))
  return result

def biasVMassCompareSigmas(data,data2,useRefFuncs,outputDir):

  refPdfNameList = data['meta']['refPdfNameList']
  pdfAltNamesDict = data['meta']['pdfAltNamesDict']
  nToys = data['meta']['nToys']
  nData = data['meta']['nData']
  catName = data['meta']['catName']
  energyStr = data['meta']['energyStr']
  sigMasses = data['meta']['sigMasses']
  sigInject = data['meta']['sigInjectNsigma']

  assert(refPdfNameList == data2['meta']['refPdfNameList'])
  assert(pdfAltNamesDict == data2['meta']['pdfAltNamesDict'])
  assert(catName == data2['meta']['catName'])
  assert(energyStr == data2['meta']['energyStr'])
  assert(sigMasses == data2['meta']['sigMasses'])
  
  if sigMasses.count(115):
    sigMasses.remove(115)
  if sigMasses.count(155):
    sigMasses.remove(155)

  ## Hack for only one alt!!
  altName = pdfAltNamesDict[refPdfNameList[0]][0]
  altTitle = PDFTITLEMAP[altName]
  gList = []
  gList2 = []
  maxNUncList = []

  canvas = root.TCanvas()

  leg = root.TLegend(0.55,0.60,0.9,0.9)
  leg.SetFillColor(0)
  leg.SetLineColor(0)
  leg.SetHeader("Reference Functions")

  maxN = 0.
  minN = 0.

  colors = [root.kBlack,root.kBlue,root.kRed,root.kGreen,root.kCyan,root.kMagenta,root.kOrange-3,root.kViolet-6]

  for iRef,refName in enumerate(refPdfNameList):
    refTitle = PDFTITLEMAP[refName]
    g = root.TGraph()
    g2 = root.TGraph()
    g.SetMarkerColor(colors[iRef])
    g.SetLineColor(colors[iRef])
    g2.SetMarkerColor(colors[iRef])
    g2.SetLineColor(colors[iRef])
    g2.SetLineStyle(2)
    maxNUnc = 0.
    for iPoint,hmass in enumerate(sigMasses):
      #print refName
      #print hmass
      #print altName
      #print iPoint
      ## Use H->gg inspired bias measure
      ##  it absorbs low stats funny business :-)
      nList = median(numpy.array(data[refName][hmass][altName]['n'])-numpy.array(data[refName][hmass]['nTrue']))
      nList2 = median(numpy.array(data2[refName][hmass][altName]['n'])-numpy.array(data2[refName][hmass]['nTrue']))
      n = median(nList)
      n2 = median(nList2)
      g.SetPoint(iPoint,hmass,n)
      g2.SetPoint(iPoint,hmass,n2)
      maxN = max(maxN,n)
      minN = min(minN,n)
      maxN = max(maxN,n2)
      minN = min(minN,n2)
      maxNUnc = max(maxNUnc,abs(n))
    gList.append(g)
    gList2.append(g2)
    leg.AddEntry(g,refTitle,"lp")
    maxNUncList.append(maxNUnc)

  gZero = root.TGraph()
  gZero.SetLineStyle(2)
  gZero.SetPoint(0,115,0.)
  gZero.SetPoint(1,155,0.)

  gError = root.TGraphErrors()
  gError.SetLineColor(0)
  gError.SetFillColor(root.kGray)
  for iPoint,hmass in enumerate(sigMasses):
    n = ERRORSDICT[energyStr][catName][hmass]
    gError.SetPoint(iPoint,hmass,0.)
    gError.SetPointError(iPoint,0.,n)
    maxN = max(maxN,n)
    minN = min(minN,-n)
  leg.AddEntry(gError,"Stat. Error","f")

  maxN *= 3.5
  minN *= 1.5
  axisHist = root.TH2F("axisHist","",1,115,155,1,minN,maxN)
  setHistTitles(axisHist,"m_{H} [GeV/c^{2}]","N_{bias}")
  axisHist.Draw()
  gError.Draw('3')
  gZero.Draw("L")
  for g in reversed(gList):
    g.Draw("PL")
  for g2 in reversed(gList2):
    g2.Draw("L")

  leg.Draw()

  energyTitle = energyStr
  energyTitle = "#sqrt{s} = "+energyTitle.replace("TeV"," TeV")
  lumiTitle = "L = {0:0.1f} fb^{{-1}}".format(lumiDict[energyStr])
  tlatex = drawStandardCaptions(canvas,TITLEMAP[catName],energyTitle,lumiTitle,preliminaryString="CMS Internal")
  tlatex.SetTextColor(root.kRed)
  tlatex.SetTextFont(62)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.76,"Solid: 0#sigma")
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.70,"Dashed: 3#sigma")

  canvas.RedrawAxis()
  saveAs(canvas,outputDir+"BiasVMassSigCompare_"+catName+"_"+energyStr)

if __name__ == "__main__":

  outputDir = "output/"
  
  for energyStr in ["7TeV","8TeV"]:
    all0SigData = getBiasStudies("output/0sig/*"+energyStr+"*.pkl")
    all3SigData = getBiasStudies("output/3sig/*"+energyStr+"*.pkl")
    for i,j in zip(all0SigData,all3SigData):
      assert(i['meta']['catName']==j['meta']['catName'])
      biasVMassCompareSigmas(i,j,None,outputDir)
