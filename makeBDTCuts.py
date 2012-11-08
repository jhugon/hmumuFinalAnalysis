#!/usr/bin/env python

from xsec import *
from helpers import *
import ROOT as root
import os
import os.path
import sys
import glob

import makeCards

def getSigLimits(fn, histName, quantile=0.99):
  f = root.TFile(fn)
  h = f.Get(histName)
  l = getMedianAndQuantileInterval(h,quantile)
  f.Close()
  return l[2],l[0]

def makeSillyHist(inHist,outHist,cutVal):
      cutBin = inHist.GetYaxis().FindBin(cutVal)
      nBinsY = inHist.GetYaxis().GetNbins()
      nBinsX = inHist.GetXaxis().GetNbins()
      #print("max bin: {0}, cutBin: {1}".format(nBinsY,cutBin))
      for iX in range(0,nBinsX+2):
        mySum = 0.
        for iY in range(cutBin,nBinsY+2):
          mySum += inHist.GetBinContent(iX,iY)
        outHist.SetBinContent(iX,mySum)

def makeSillyFiles(inDir,outDir, runperiod,
            rebinY=100
    ):
  incLimits = getSigLimits(dataDir+incSigFilename,"BDTHistMuonOnly")
  vbfLimits = getSigLimits(dataDir+vbfSigFilename,"BDTHistVBF")
  print("Inc Limits: {0} \tVBF Limits: {1}".format(incLimits,vbfLimits))
  incHistname = "BDTHistMuonOnlyVMass"
  vbfHistname = "BDTHistVBFVMass"
  filenames = glob.glob(dataDir+"/*_"+RUNPERIOD+".root")
  for fn in filenames:
    infile = root.TFile(fn)
    inc = infile.Get(incHistname)
    vbf = infile.Get(vbfHistname)
    inc.Rebin2D(1,rebinY)
    vbf.Rebin2D(1,rebinY)
    out = infile.Get("mDiMu")
  
    #outFileName = os.path.splitext(fn)[0]
    outFileName = os.path.basename(fn)
    outFileName = outDir +"/"+ outFileName
    outFileName = os.path.normpath(outFileName)
    print outFileName
    outFile = root.TFile(outFileName,"RECREATE")
    outFile.cd()
  
    for iCut in range(0,10):
      incCutVal = incLimits[0]+iCut*(incLimits[1]-incLimits[0])/10.
      vbfCutVal = vbfLimits[0]+iCut*(vbfLimits[1]-vbfLimits[0])/10.
      outVBF = out.Clone("mDiMuBDTVBFCut{0:.2f}".format(incCutVal))
      outInc = out.Clone("mDiMuBDTIncCut{0:.2f}".format(vbfCutVal))
      outVBF.Reset()
      outInc.Reset()
      makeSillyHist(inc,outInc,incCutVal)
      makeSillyHist(vbf,outVBF,vbfCutVal)
      print("Inc cut: {0:.2f}, sum: {1} \tVBF cut: {2:.2f}, sum: {3}".format(incCutVal,outInc.Integral(),vbfCutVal,outVBF.Integral()))
      outInc.Write()
      outVBF.Write()
    outFile.Close()

if __name__ == "__main__":

  RUNPERIOD = "7TeV"

  dataDir = "input/"
  outDir = "bdtCutInputs/"
  incSigFilename = "ggHmumu125_"+RUNPERIOD+".root"
  vbfSigFilename = "vbfHmumu125_"+RUNPERIOD+".root"

  makeSillyFiles(dataDir,outDir,RUNPERIOD)
