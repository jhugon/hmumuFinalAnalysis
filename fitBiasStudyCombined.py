#! /usr/bin/env python

import datetime
import sys
import os
import os.path
import re
import math
import cPickle
import glob
import ROOT as root
root.gROOT.SetBatch(True)

from helpers import *
from xsec import *

import scipy.stats
from numpy import mean, median, corrcoef, percentile,vstack
from numpy import std as stddev

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT


################################################################################################
################################################################################################
################################################################################################

def setYMaxAndDrawVertLines(hist,x=None):
  ymax = 0.
  for i in range(1,hist.GetXaxis().GetNbins()+1):
    ymax = max(ymax,hist.GetBinContent(i))
  ymax * 1.5
  hist.GetYaxis().SetRangeUser(0.,ymax*1.5)
  if x == None:
      return
  line = root.TLine()
  line.SetLineColor(root.kRed)
  line.SetLineWidth(3)
  line.DrawLine(x,0,x,ymax*1.05)
  return line

def getData(refDir,altDir):
  refFns = glob.glob(refDir+"*.muToys.root")
  altFns = glob.glob(altDir+"*.muToys.root")
  refFns = sorted(refFns)
  altFns = sorted(altFns)
  refJobNs = []
  prefix="CombSplitAll_7P8TeV_125.0.txt"
  for refFn in refFns:
    refMatch = re.search(prefix+r"\.([0-9]+)\.muToys\.root",refFn)
    if refMatch:
      refJobNs.append(refMatch.group(1))
    else:
      refJobNs.append(-1)
  altJobNs = []
  for altFn in altFns:
    altMatch = re.search(prefix+r"\.([0-9]+)\.muToys\.root",altFn)
    if altMatch:
      altJobNs.append(altMatch.group(1))
    else:
      altJobNs.append(-1)

  altJobNsSet = set(altJobNs)
  refJobNsSet = set(refJobNs)

  refDifSet = refJobNsSet.difference(altJobNsSet)
  altDifSet = altJobNsSet.difference(refJobNsSet)

  for i in reversed(range(len(altFns))):
    n = altJobNs[i]
    if n == -1 or n in altDifSet:
      altFns.pop(i)
      altJobNs.pop(i)

  for i in reversed(range(len(refFns))):
    n = refJobNs[i]
    if n == -1 or n in refDifSet:
      refFns.pop(i)
      refJobNs.pop(i)

  assert(len(altFns)==len(refFns))
  assert(len(altFns)==len(altJobNs))
  assert(len(refFns)==len(refJobNs))

  result = []

  for i in range(len(altFns)):
    altFn = altFns[i]
    refFn = refFns[i]
    altN = altJobNs[i]
    refN = refJobNs[i]
    assert(altN==refN)
    altTF = root.TFile(altFn)
    refTF = root.TFile(refFn)
    altTree = altTF.Get("limit")
    refTree = refTF.Get("limit")

    # Figure out poper iToys
    iToySetAlt = set()
    for iLimit in range(3,altTree.GetEntries(),4):
      altTree.GetEntry(iLimit)
      iToySetAlt.add(altTree.iToy)
    iToySetRef = set()
    for iLimit in range(3,refTree.GetEntries(),4):
      refTree.GetEntry(iLimit)
      iToySetRef.add(refTree.iToy)
    iToySetBoth = iToySetAlt.intersection(iToySetRef)
    for iToy in iToySetBoth:
      for iLimit in range(3,refTree.GetEntries(),4):
        refTree.GetEntry(iLimit)
        if refTree.iToy == iToy:
            break
      for iLimit in range(3,altTree.GetEntries(),4):
        altTree.GetEntry(iLimit)
        if altTree.iToy == iToy:
            break
      assert(altTree.iToy==refTree.iToy)
      assert(altTree.iSeed==refTree.iSeed)
      #if altTree.limitErr > 4. or altTree.limitErr < 1.:
      #  continue
      #if refTree.limitErr > 4. or refTree.limitErr < 1.:
      #  continue
      point = {}
      point['muRef'] = refTree.limit
      point['muAlt'] = altTree.limit
      point['muErrRef'] = refTree.limitErr
      point['muErrAlt'] = altTree.limitErr
      result.append(point)
  return result
      

if __name__ == "__main__":
  inputRefDir= "/afs/cern.ch/user/j/jhugon/work/private/stats/CMSSW_6_1_1/finalAnalysisVoigtExp110160/statsCards/"
  inputAltDir= "/afs/cern.ch/user/j/jhugon/work/private/stats/CMSSW_6_1_1/finalAnalysis/statsCards/"
  outputPrefix = "statsOutput/biasComb_"

  data = getData(inputRefDir,inputAltDir)

  refMus = [toy['muRef'] for toy in data]
  altMus = [toy['muAlt'] for toy in data]
  refErrMus = [toy['muErrRef'] for toy in data]
  altErrMus = [toy['muErrAlt'] for toy in data]
  refZs = [toy['muRef']/toy['muErrRef'] for toy in data]
  altZs = [toy['muAlt']/toy['muErrAlt'] for toy in data]

  pulls = [(toy['muAlt']-toy['muRef'])/toy['muErrRef'] for toy in data]

  PRELIMINARYSTRING="CMS Internal"
  
  canvas = root.TCanvas()
  tlatex = root.TLatex()
  tlatex.SetNDC()
  tlatex.SetTextFont(root.gStyle.GetLabelFont())
  tlatex.SetTextSize(0.04)
  caption = "H#rightarrow#mu#mu Combination"
  caption2 = ""
  caption3 = ""
  caption4 = ""
  iHist = 0

  hmass = 125.
  refPdfName="Old"
  pdfAltName="ExpMOverSq"

  ## Pull Plot
  hist = root.TH1F("hist"+str(iHist),"",60,-3,3)
  setHistTitles(hist,"(N_{sig}(Alt)-N_{sig}(Ref))/#DeltaN_{sig}(Alt)","N_{Toys}")
  iHist += 1
  for pull in pulls:
      hist.Fill(pull)
  hist.Draw()
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfName])
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+PDFTITLEMAP[pdfAltName])
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.60,"Median: {0:.2f}".format(median(pulls)))
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
  line = setYMaxAndDrawVertLines(hist,median(pulls))
  canvas.RedrawAxis()
  saveAs(canvas,outputPrefix+"_"+str(hmass)+"_Pulls_Ref"+refPdfName+"_Alt"+pdfAltName)
  canvas.Clear()

  ## Mu Ref Plot
  hist = root.TH1F("hist"+str(iHist),"",60,-10,10)
  setHistTitles(hist,"#mu Reference","N_{Toys}")
  iHist += 1
  for mu in refMus:
      hist.Fill(mu)
  hist.Draw()
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfName])
  #tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+PDFTITLEMAP[pdfAltName])
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Median: {0:.2f}".format(median(refMus)))
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Mean: {0:.2f}".format(mean(refMus)))
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"#sigma: {0:.2f}".format(stddev(refMus)))
  line = setYMaxAndDrawVertLines(hist,None)
  canvas.RedrawAxis()
  saveAs(canvas,outputPrefix+"_"+str(hmass)+"_MuRef_Ref"+refPdfName)
  canvas.Clear()

  ## Mu Alt Plot
  hist = root.TH1F("hist"+str(iHist),"",60,-10,10)
  setHistTitles(hist,"#mu Alternate","N_{Toys}")
  iHist += 1
  for mu in altMus:
      hist.Fill(mu)
  hist.Draw()
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfName])
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+PDFTITLEMAP[pdfAltName])
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Median: {0:.2f}".format(median(altMus)))
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Mean: {0:.2f}".format(mean(altMus)))
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"#sigma: {0:.2f}".format(stddev(altMus)))
  line = setYMaxAndDrawVertLines(hist,None)
  canvas.RedrawAxis()
  saveAs(canvas,outputPrefix+"_"+str(hmass)+"_MuAlt_Ref"+refPdfName+"_Alt"+pdfAltName)
  canvas.Clear()

  ## MuErr Ref Plot
  hist = root.TH1F("hist"+str(iHist),"",30,0,10)
  setHistTitles(hist,"#Delta#mu Reference","N_{Toys}")
  iHist += 1
  for mu in refErrMus:
      hist.Fill(mu)
  hist.Draw()
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfName])
  #tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+PDFTITLEMAP[pdfAltName])
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Median: {0:.2f}".format(median(refErrMus)))
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Mean: {0:.2f}".format(mean(refErrMus)))
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"#sigma: {0:.2f}".format(stddev(refErrMus)))
  line = setYMaxAndDrawVertLines(hist,None)
  canvas.RedrawAxis()
  saveAs(canvas,outputPrefix+"_"+str(hmass)+"_MuErrRef_Ref"+refPdfName)
  canvas.Clear()

  ## MuErr Alt Plot
  hist = root.TH1F("hist"+str(iHist),"",30,0,10)
  setHistTitles(hist,"#Delta#mu Alternate","N_{Toys}")
  iHist += 1
  for mu in altErrMus:
      hist.Fill(mu)
  hist.Draw()
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfName])
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+PDFTITLEMAP[pdfAltName])
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Median: {0:.2f}".format(median(altErrMus)))
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Mean: {0:.2f}".format(mean(altErrMus)))
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"#sigma: {0:.2f}".format(stddev(altErrMus)))
  line = setYMaxAndDrawVertLines(hist,None)
  canvas.RedrawAxis()
  saveAs(canvas,outputPrefix+"_"+str(hmass)+"_MuErrAlt_Ref"+refPdfName+"_Alt"+pdfAltName)
  canvas.Clear()

  ## Z Ref Plot
  hist = root.TH1F("hist"+str(iHist),"",60,-10,10)
  setHistTitles(hist,"#mu/#Delta#mu Reference","N_{Toys}")
  iHist += 1
  for mu in refZs:
      hist.Fill(mu)
  hist.Draw()
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfName])
  #tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+PDFTITLEMAP[pdfAltName])
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Median: {0:.2f}".format(median(refZs)))
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Mean: {0:.2f}".format(mean(refZs)))
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"#sigma: {0:.2f}".format(stddev(refZs)))
  line = setYMaxAndDrawVertLines(hist,None)
  canvas.RedrawAxis()
  saveAs(canvas,outputPrefix+"_"+str(hmass)+"_Z_Ref"+refPdfName)
  canvas.Clear()

  ## Z Alt Plot
  hist = root.TH1F("hist"+str(iHist),"",60,-10,10)
  setHistTitles(hist,"#mu/#Delta#mu Alternate","N_{Toys}")
  iHist += 1
  for mu in altZs:
      hist.Fill(mu)
  hist.Draw()
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfName])
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+PDFTITLEMAP[pdfAltName])
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Median: {0:.2f}".format(median(altZs)))
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Mean: {0:.2f}".format(mean(altZs)))
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"#sigma: {0:.2f}".format(stddev(altZs)))
  line = setYMaxAndDrawVertLines(hist,None)
  canvas.RedrawAxis()
  saveAs(canvas,outputPrefix+"_"+str(hmass)+"_Z_Ref"+refPdfName+"_Alt"+pdfAltName)
  canvas.Clear()

  ## Mu Correlation Plot
  hist = root.TH2F("hist"+str(iHist),"",30,-10,10,30,-10,10)
  setHistTitles(hist,"#mu Alternate","#mu Reference")
  iHist += 1
  #for muA,muR in zip(altMus,refMus):
  #    hist.Fill(muA,muR)
  #hist.Draw("col")
  graph = root.TGraph()
  for iPoint,muA,muR in zip(range(len(altMus)),altMus,refMus):
      graph.SetPoint(iPoint,muA,muR)
  hist.Draw("")
  graph.Draw("P")
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfName])
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+PDFTITLEMAP[pdfAltName])
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.60,"#rho =  {0:.3f}".format(corrcoef(altMus,refMus)[0,1]))
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
  #tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Median: {0:.2f}".format(median(altZs)))
  #tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Mean: {0:.2f}".format(mean(altZs)))
  #tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"#sigma: {0:.2f}".format(stddev(altZs)))
  #line = setYMaxAndDrawVertLines(hist,None)
  canvas.RedrawAxis()
  saveAs(canvas,outputPrefix+"_"+str(hmass)+"_Correlation_Ref"+refPdfName+"_Alt"+pdfAltName)
  canvas.Clear()


  print "\n#########################################\n"
  print "NToys: ",len(data)
  print "\n#########################################\n"
  print "Reference Median Mu: ",median(refMus)
  print "Reference Mean Mu: ",mean(refMus)
  print "Reference Sigma Mu: ",stddev(refMus)
  print "Alternate Median Mu: ",median(altMus)
  print "Alternate Mean Mu: ",mean(altMus)
  print "Alternate Sigma Mu: ",stddev(altMus)
  print "\n#########################################\n"
  print "Reference Median ErrMu: ",median(refErrMus)
  print "Reference Mean ErrMu: ",mean(refErrMus)
  print "Reference Sigma ErrMu: ",stddev(refErrMus)
  print "Alternate Median ErrMu: ",median(altErrMus)
  print "Alternate Mean ErrMu: ",mean(altErrMus)
  print "Alternate Sigma ErrMu: ",stddev(altErrMus)
  print "\n#########################################\n"
  print "Mu Alt & Mu Ref Correlation:",corrcoef(altMus,refMus)[0,1]
  print "\n#########################################\n"
  print "Reference Median Z: ",median(refZs)
  print "Reference Mean Z: ",mean(refZs)
  print "Reference Sigma Z: ",stddev(refZs)
  print "Alternate Median Z: ",median(altZs)
  print "Alternate Mean Z: ",mean(altZs)
  print "Alternate Sigma Z: ",stddev(altZs)
  print "\n#########################################\n"
  print "Median Bias: ",median(pulls)
  print "\n#########################################\n"
  
