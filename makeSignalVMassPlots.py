#!/usr/bin/env python

import ROOT as root
from helpers import *
from xsec import *
import datetime
import sys
import os.path
import glob
import makeCards
from makeCards import makePDFSig, SIGNALFIT

####################################################

class SigVMass:
  def __init__(self,infilenames,categories,histSuffix="mDiMu"):
    self.infilenames = infilenames
    self.categories = categories
    self.histSuffix = histSuffix
    self.files = []
    self.masses = []
    for i in infilenames:
      self.files.append(root.TFile(i))
      mass = -1.
      match = re.match(r".*Hmumu([\d.]+)_.*",i)
      if match:
        mass = float(match.group(1))
      self.masses.append(mass)
    self.hists = {}
    for c in categories:
      tmpList = []
      for i in self.files:
        tmpList.append(i.Get(os.path.normpath(c+"/"+histSuffix)))
      self.hists[c] = tmpList

    self.workspace = root.RooWorkspace("w")
    wImport = getattr(self.workspace,"import")
    self.wImport = wImport

    self.params = {}
    for c in categories:
      tmpList = []
      for inHist,mass in zip(self.hists[c],self.masses):
        inHist.Rebin(1)
        mMuMu = root.RooRealVar("mMuMu","mMuMu",mass-20.0,mass+20.0)
        mMuMu.setRange("signalfit",mass-15.0,mass+15.0)
        #makePDFSig = makeCards.makePDFSigGaus
        #makePDFSig = makeCards.makePDFSigCBConvGaus
        params,debug = makePDFSig(c+str(mass),inHist,mMuMu,
                            SIGNALFIT[0],SIGNALFIT[1],
                            wImport,c+str(mass)
                            )
        tmpList.append(params)
      self.params[c] = tmpList

    self.eff = {}
    tEff = root.TEfficiency()
    for c in categories:
      tmpList = []
      tmpListErrUp = []
      tmpListErrDown = []
      for inHist,mass,inFn in zip(self.hists[c],self.masses,self.infilenames):
        tmpSampleName = os.path.splitext(os.path.split(inFn)[1])[0]
        print tmpSampleName
        tmpEff = inHist.Integral()/nEventsMap[tmpSampleName]
        cLevel = 0.683
        tmpList.append((
            tEff.ClopperPearson(nEventsMap[tmpSampleName],int(inHist.Integral()),cLevel,False),
            tmpEff,
            tEff.ClopperPearson(nEventsMap[tmpSampleName],int(inHist.Integral()),cLevel,True)
            )
            )
      self.eff[c] = tmpList

  def __str__(self):
    result = "============================================\n"
    result += "SigVMass\n"
    result += "============================================\n"
    result += "Hist Suffix: "+self.histSuffix+"\n"
    result += "Categories:\n"
    for c in self.categories:
      result += "  "+c +"\n"
    result += "\n{0:<6} {1:<}\n".format("Mass","Filename")
    result += "-------------------------------------------\n"
    for m,i in zip(self.masses,self.infilenames):
      result += "{0:<6} {1:<}\n".format(m,i)
    result += "\n"
    if len(self.infilenames) == 0:
        return result
    for c in self.categories:
      result += c +"\n"
      result += "-------------------------------------------\n"
      for iparam in range(len(self.params[c][0])):
        result += re.match(r".*_([a-zA-Z0-9]+)$",self.params[c][0][iparam].name).group(1) + "\n"
        for iFile, mass in zip(range(len(self.files)),self.masses):
          tmpParam = self.params[c][iFile][iparam]
          result += "  {0:<6} {1:<10.5g} +/- {2:<10.5g}\n".format(mass,tmpParam.nominal,
                            max(abs(tmpParam.highErr),abs(tmpParam.lowErr)))

      result += "\nEfficiency: \n"
      for iFile, mass in zip(range(len(self.files)),self.masses):
          tmpEff = self.eff[c][iFile]
          result += "  {0:<6} {1:<10.3%} + {2:<10.3%} - {3:<10.3%}\n".format(mass,tmpEff[1],
                            tmpEff[2]-tmpEff[1],tmpEff[1]-tmpEff[0])
    return result + "\n"

  def plots(self,outDir):
    if len(self.infilenames) == 0:
        return
    canvas = root.TCanvas("canvas")
    graphs = []
    for c in self.categories:
      for iparam in range(len(self.params[c][0])):
        canvas.Clear()
        graph = root.TGraphAsymmErrors()
        paramName =  re.match(r".*_([a-zA-Z0-9]+)$",self.params[c][0][iparam].name).group(1)
        for iFile, mass in zip(range(len(self.files)),self.masses):
          tmpParam = self.params[c][iFile][iparam]
          graph.SetPoint(iFile,mass,tmpParam.nominal)
          graph.SetPointError(iFile,0.0,0.0,tmpParam.lowErr,tmpParam.highErr)
        setHistTitles(graph,"Generator m_{H} [GeV/c^{2}]",paramName)
        if paramName == "mix":
          graph.GetYaxis().SetRangeUser(0,1)
        if paramName == "Width" or paramName == "Width2":
          graph.GetYaxis().SetRangeUser(0,6)
        graph.Draw("ap")
        saveAs(canvas,outDir+"_"+c+"_"+paramName)
      # Efficiency
      canvas.Clear()
      graph = root.TGraphAsymmErrors()
      for iFile, mass in zip(range(len(self.files)),self.masses):
          tmpEff = self.eff[c][iFile]
          graph.SetPoint(iFile,mass,tmpEff[1])
          graph.SetPointError(iFile,0.0,0.0,tmpEff[1]-tmpEff[0],tmpEff[2]-tmpEff[1])
      setHistTitles(graph,"Generator m_{H} [GeV/c^{2}]","Acceptance #times Efficiency")
      graph.GetYaxis().SetRangeUser(self.eff[c][iFile][0]*0.9,self.eff[c][iFile][2]*1.1)
      graph.Draw("ap")
      saveAs(canvas,outDir+"_"+c+"_Eff")
  
if __name__ == "__main__":
  root.gErrorIgnoreLevel = root.kWarning
  root.gROOT.SetBatch(True)
  setStyle()

  indir = "input/preApproveSample/"
  outdir = "output/"

  categories = ["IncPreselPtG10","VBFBDTCut"]

  for e in ["8TeV"]:
    for p in ["vbf","gg"]:
      filenames = glob.glob(indir+p+"Hmumu*"+e+".root")
      svm = SigVMass(filenames,categories)
      svm.plots(outdir+p+e)
      print svm
