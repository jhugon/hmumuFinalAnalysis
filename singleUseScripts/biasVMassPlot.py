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

def getBiasStudies():
  """
  loads in biasstudy raw .pkl files
  returns a list of combined dictionaries
  """
  result = []
  inputPklFiles = glob.glob("output/"+"*.pkl")
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

def printBiasVMass(data,useRefFuncs,outputDir):

  refPdfNameList = data['meta']['refPdfNameList']
  pdfAltNamesDict = data['meta']['pdfAltNamesDict']
  nToys = data['meta']['nToys']
  nData = data['meta']['nData']
  catName = data['meta']['catName']
  energyStr = data['meta']['energyStr']
  sigMasses = data['meta']['sigMasses']
  sigInject = data['meta']['sigInjectNsigma']

  if sigMasses.count(115):
    sigMasses.remove(115)
  if sigMasses.count(155):
    sigMasses.remove(155)

  ## Hack for only one alt!!
  altName = pdfAltNamesDict[refPdfNameList[0]][0]
  altTitle = PDFTITLEMAP[altName]
  gList = []
  maxNUncList = []

  canvas = root.TCanvas()

  leg = root.TLegend(0.58,0.65,0.9,0.9)
  leg.SetFillColor(0)
  leg.SetLineColor(0)
  leg.SetHeader("Reference Functions")

  maxN = 0.
  minN = 0.

  colors = [root.kBlack,root.kBlue,root.kRed,root.kGreen,root.kCyan,root.kMagenta,root.kOrange-3,root.kViolet-6]

  for iRef,refName in enumerate(refPdfNameList):
    refTitle = PDFTITLEMAP[refName]
    g = root.TGraph()
    g.SetMarkerColor(colors[iRef])
    g.SetLineColor(colors[iRef])
    maxNUnc = 0.
    for iPoint,hmass in enumerate(sigMasses):
      #print refName
      #print hmass
      #print altName
      #print iPoint
      ## Use H->gg inspired bias measure
      ##  it absorbs low stats funny business :-)
      nList = median(numpy.array(data[refName][hmass][altName]['n'])-numpy.array(data[refName][hmass]['nTrue']))
      n = median(nList)
      g.SetPoint(iPoint,hmass,n)
      maxN = max(maxN,n)
      minN = min(minN,n)
      maxNUnc = max(maxNUnc,abs(n))
    gList.append(g)
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

  leg.Draw()

  energyTitle = energyStr
  energyTitle = "#sqrt{s} = "+energyTitle.replace("TeV"," TeV")
  tlatex = drawStandardCaptions(canvas,"Bias for "+altTitle,TITLEMAP[catName],energyTitle,preliminaryString="CMS Internal")

  canvas.RedrawAxis()
  saveAs(canvas,"BiasVMass_"+catName+"_"+energyStr+"_sig"+str(int(sigInject)))

  ### Now to print table
  statUnc125 = ERRORSDICT[energyStr][catName][125]
  sigInjectTitle = r"${0}\sigma$ Signal Injected".format(sigInject)
  if sigInject == 0.:
    sigInjectTitle = "No Signal Injected"
  outStr = r"\normalsize"+"\n"
  outStr += r"{0} {1} {2} \\".format(TITLEMAP[catName],energyStr.replace("TeV"," TeV"), sigInjectTitle)+"\n"
  outStr += "Bias for "+altName+r" \\ "+"\n"
  outStr += r"\tiny"+"\n"
  outStr += r"\begin{tabular}{|l|r|r|} \hline"+"\n"
  outStr += r"Reference PDF & N_{bias} & Relative Bias \\ \hline \hline"+"\n"
  for refPdfName,maxNUnc in zip(refPdfNameList,maxNUncList):
    refPdfTitle = PDFTITLEMAP[refPdfName]
    outStr += r"{0:40} & {1:15.1f} & {2:15.0f}\% \\ \hline".format(refPdfTitle,maxNUnc,maxNUnc*100./statUnc125)+"\n"
  outStr += r"\end{tabular}"+"\n\n"
  print outStr

  resultNBiasDict = {}
  for refPdfName,maxNUnc in zip(refPdfNameList,maxNUncList):
    resultNBiasDict[refPdfName] = maxNUnc
  return resultNBiasDict

if __name__ == "__main__":
  onlyVoitRefs = ["Old","VoigtPMm2","VoigtPExpMm2"]
  voitAndSMRefs = ["Old","ExpMOverSq","VoigtPMm2","VoigtPExpMm2"]
  allRefs = ["Old","ExpMOverSq","SumExp","VoigtPMm2","Bernstein","VoigtPExpMm2"]
  bernRefs = ["3Bernstein","4Bernstein","5Bernstein"]

  outputDir = "output/"

  nBiasMaxDict = {}
  energyStr = None
  sigInject = None
  
  allData = getBiasStudies()
  for i in allData:
    nBiasMaxDict[i['meta']['catName']] = printBiasVMass(i,None,outputDir)
    energyStr = i['meta']['energyStr']
    sigInject = i['meta']['sigInjectNsigma']
  assert(energyStr)
  if sigInject == None:
    sigInject = 0.
  outPklFile = open("BiasDataVRefFunc_{0}_{1:.0f}Sig.pkl".format(energyStr,sigInject),"w")
  cPickle.dump(nBiasMaxDict,outPklFile)
  outPklFile.close()
