#! /usr/bin/env python

from ROOT import gSystem

import singleHelpers
import datetime
import sys
import os
import re
import math
import cPickle
import ROOT as root
root.gSystem.Load('libRooFit')
root.gROOT.SetBatch(True)
import scipy.stats

from helpers import *
from makeCards import *
from xsec import *
from fitOrderChooser import makePDFBakSumExp

from numpy import mean, median, corrcoef, percentile
from numpy import std as stddev

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

nToys = 500
nBins1d = 50
nBins2d = 20

inputFile = root.TFile("output/debug_RooFit_Bernstein_Jets01PassPtG10BB_8TeV_job0.root")

w = inputFile.Get("debugWorkspace")
#w.Print()

dimuonMass = w.var("dimuonMass")
nSig = w.var("nSig")

## Load PDF from file
pdfs =  rooArgSet2List(w.allPdfs())
mssmPdf = None
bernPdf = None
sigPdf = None
for i in pdfs:
  if "MSSM" in i.GetName():
    mssmPdf = i
  if "Bernstein" == i.GetName():
    bernPdf = i
  if "sigPDF" in i.GetName() and i.GetName()[-1] == "E":
    sigPdf = i

realData = w.data("realDataJets01PassPtG10BB8TeV")

frRealData = bernPdf.fitTo(realData,
                 root.RooFit.Minos(True),
                 PRINTLEVEL,
                 root.RooFit.Save()
    )

nBak = root.RooRealVar("nBak","N_{bkg}",realData.sumEntries(),nSig.getMax()*2,realData.sumEntries()*2)
mssmPdfE = root.RooExtendPdf("mssmPdfE","Extended Background PDF",mssmPdf,nBak)
mssmPdfSB = root.RooAddPdf("mssmPdfSB","S+B PDF",root.RooArgList(mssmPdfE,sigPdf))

bernPdfE = root.RooExtendPdf("bernPdfE","Extended Background PDF",bernPdf,nBak)
bernPdfSB = root.RooAddPdf("bernPdfSB","S+B PDF",root.RooArgList(bernPdfE,sigPdf))


params = rooArgSet2List(mssmPdfSB.getParameters(realData))
fitVals = {}
hists1D = {}
for param in params:
  if not param.isConstant():
    name = param.GetName()
    shortName = name
    if '_' in name:
      shortName = name[-[i for i in reversed(name)].index("_"):]
    fitVals[name] = []
    pMin = param.getMin()
    pMax = param.getMax()
    if shortName == "nSig":
      pMin = -250
      pMax = 250
    hists1D[name] = root.TH1F("1d_"+shortName,"",nBins1d,pMin,pMax)
    hists1D[name].Sumw2()
    setHistTitles(hists1D[name],shortName,"N_{toys}")

nSigRefVals = []

pNames = sorted(fitVals.keys())
canvas = root.TCanvas("canvas")

for iToy in range(nToys):
  toyData = bernPdf.generateBinned(root.RooArgSet(dimuonMass),realData.sumEntries())
  fr = mssmPdfSB.fitTo(toyData,
                 root.RooFit.Minos(True),
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  floatPars = rooArgSet2Dict(fr.floatParsFinal())
  for pName in pNames:
    fitVals[pName].append( floatPars[pName].getVal() )

  fr = bernPdfSB.fitTo(toyData,
                 root.RooFit.Minos(True),
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  floatPars = rooArgSet2Dict(fr.floatParsFinal())
  nSigRefVals.append(floatPars["nSig"].getVal())

tlatex = root.TLatex()
tlatex.SetNDC()
tlatex.SetTextFont(root.gStyle.GetLabelFont())
tlatex.SetTextSize(0.04)
  
for pName in pNames:
  for val in fitVals[pName]:
    hists1D[pName].Fill(val)

  hists1D[pName].Draw()
  if pName == "nSig":
    fitFn = root.TF1("gausForFit"+pName,"gaus",hists1D[pName].GetXaxis().GetBinLowEdge(1),hists1D[pName].GetXaxis().GetBinUpEdge(nBins1d))
    fitResult = hists1D[pName].Fit(fitFn,"LSMEQ") 
    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Mean: {0:.2f} #pm {1:.2f}".format(fitFn.GetParameter(1),fitFn.GetParError(1)))
    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.79,"#sigma: {0:.2f} #pm {1:.2f}".format(fitFn.GetParameter(2),fitFn.GetParError(2)))
  name = hists1D[pName].GetName()
  drawStandardCaptions(canvas,"0,1-Jet Tight BB #sqrt{s} = 8 TeV","Toys Generated From 4-Bernstin","Fit with MSSM","nToys: "+str(nToys),preliminaryString="")
  saveAs(canvas,"output/fitInvestigationsEnsemble_"+name)

hists2D = {}
graphs2D = {}
for i in range(len(pNames)):
  iName = pNames[i]
  iNameShort = iName
  if '_' in iName:
    iNameShort = iName[-[k for k in reversed(iName)].index("_"):]
  hists2D[iName] = {}
  graphs2D[iName] = {}
  for j in range(i+1,len(pNames)):
    jName = pNames[j]
    jNameShort = jName
    if '_' in jName:
      jNameShort = jName[-[k for k in reversed(jName)].index("_"):]
    name = "2d_"+iNameShort+"_"+jNameShort
    hists2D[iName][jName] = root.TH2F(name,"",
                                            nBins2d,
                                            hists1D[iName].GetXaxis().GetBinLowEdge(1),
                                            hists1D[iName].GetXaxis().GetBinUpEdge(nBins1d),
                                            nBins2d,
                                            hists1D[jName].GetXaxis().GetBinLowEdge(1),
                                            hists1D[jName].GetXaxis().GetBinUpEdge(nBins1d),
                                        )
    setHistTitles(hists2D[iName][jName],hists1D[iName].GetXaxis().GetTitle(),hists1D[jName].GetXaxis().GetTitle())
    graphs2D[iName][jName] = root.TGraph()
    for iToy in range(nToys):
      graphs2D[iName][jName].SetPoint(iToy,fitVals[iName][iToy],fitVals[jName][iToy])
      #hists2D[iName][jName].Fill(fitVals[iName][iToy],fitVals[jName][iToy])
    hists2D[iName][jName].Draw()
    graphs2D[iName][jName].Draw("P")
    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"#rho: {0:.2f}".format(scipy.stats.pearsonr(fitVals[iName],fitVals[jName])[0]))
    drawStandardCaptions(canvas,"0,1-Jet Tight BB #sqrt{s} = 8 TeV","Toys Generated From 4-Bernstin","Fit with MSSM","nToys: "+str(nToys),preliminaryString="")
    saveAs(canvas,"output/fitInvestigationsEnsemble_"+name)
    
outpklfile = open("output/fitInvestigationsEnsemble.pkl",'w')
cPickle.dump(fitVals,outpklfile)
outpklfile.close()
    
nSigRefHist = root.TH1F("nSigRefHist","",nBins1d,-250,250)
nSigRefHist.Sumw2()
setHistTitles(nSigRefHist,"N_{sig}(Ref)","N_{toys}")
for val in nSigRefVals:
  nSigRefHist.Fill(val)
nSigRefHist.Draw()
drawStandardCaptions(canvas,"0,1-Jet Tight BB #sqrt{s} = 8 TeV","Toys Generated From 4-Bernstin","Fit with MSSM","nToys: "+str(nToys),preliminaryString="")
saveAs(canvas,"output/fitInvestigationsEnsemble_1d_nSigRef")

nBiasVals = [alt-ref for alt,ref in zip(fitVals['nSig'],nSigRefVals)]
biasHist = root.TH1F("biasHist","",nBins1d,-250,250)
biasHist.Sumw2()
setHistTitles(biasHist,"N_{sig}(Alt) - N_{sig}(Ref)","N_{toys}")
for val in nBiasVals:
  biasHist.Fill(val)
biasHist.Draw()
drawStandardCaptions(canvas,"0,1-Jet Tight BB #sqrt{s} = 8 TeV","Toys Generated From 4-Bernstin","Fit with MSSM","nToys: "+str(nToys),preliminaryString="")
saveAs(canvas,"output/fitInvestigationsEnsemble_1d_bias")
    
for pName in pNames:
  shortName = pName
  if '_' in pName:
    shortName = pName[-[i for i in reversed(pName)].index("_"):]
  name = "2d_bias_"+shortName
  axisHist = root.TH2F(name,"",
                                          nBins2d,
                                          hists1D[pName].GetXaxis().GetBinLowEdge(1),
                                          hists1D[pName].GetXaxis().GetBinUpEdge(nBins1d),
                                          nBins2d,
                                          -250,250
                                      )
  setHistTitles(axisHist,shortName,"N_{sig}(Alt)-N_{sig}(Ref)")
  graph = root.TGraph()
  for iToy,val in enumerate(fitVals[pName]):
    nSigRef = nSigRefVals[iToy]
    nBias = nBiasVals[iToy]
    graph.SetPoint(iToy,val,nBias)
  axisHist.Draw()
  graph.Draw('P')
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"#rho: {0:.2f}".format(scipy.stats.pearsonr(fitVals[pName],nBiasVals)[0]))
  drawStandardCaptions(canvas,"0,1-Jet Tight BB #sqrt{s} = 8 TeV","Toys Generated From 4-Bernstin","Fit with MSSM","nToys: "+str(nToys),preliminaryString="")
  saveAs(canvas,"output/fitInvestigationsEnsemble_"+name)

