#!/usr/bin/env python

from singleHelpers import *
from helpers import *
import ROOT as root
import cPickle
from xsec import *
from scipy.stats import norm


PRELIMINARYSTRING="CMS Internal"
root.gStyle.SetMarkerSize(0.9)

###############################################################33

limitsBaselineF = open("limitsJets01PassPtG10BB_8TeV_Baseline.pkl")
limits20BaselineF = open("limitsJets01PassPtG10BB_8TeV_20GeVBaseline.pkl")
limits20Bern4F = open("limitsJets01PassPtG10BB_8TeV_20GeV4Bernstein.pkl")
limits20Bern2F = open("limitsJets01PassPtG10BB_8TeV_20GeV2Bernstein.pkl")

limitsBaseline = cPickle.load(limitsBaselineF)
limits20Baseline = cPickle.load(limits20BaselineF)
limits20Bern4 = cPickle.load(limits20Bern4F)
limits20Bern2 = cPickle.load(limits20Bern2F)

obsBaselineG = root.TGraph()
obs20BaselineG = root.TGraph()
obs20Bern4G = root.TGraph()
obs20Bern2G = root.TGraph()
expBaselineG = root.TGraph()
exp20BaselineG = root.TGraph()
exp20Bern4G = root.TGraph()
exp20Bern2G = root.TGraph()

oneSigBaselineG = root.TGraphAsymmErrors()
oneSigBaselineG.SetFillColor(root.kGreen)
oneSigBaselineG.SetLineStyle(0)
twoSigBaselineG = root.TGraphAsymmErrors()
twoSigBaselineG.SetFillColor(root.kYellow)
twoSigBaselineG.SetLineStyle(0)

for i,entry in zip(range(len(limitsBaseline)),limitsBaseline):

  mass = float(entry[0])

  if (mass==155):
    continue
  obsBaselineG.SetPoint(i,mass,float(entry[1]))

  median = float(entry[4])
  low2sig = median - float(entry[2])
  low1sig = median - float(entry[3])
  high1sig = float(entry[5]) - median
  high2sig = float(entry[6]) - median

  expBaselineG.SetPoint(i,mass,median)
  oneSigBaselineG.SetPoint(i,mass,median) 
  twoSigBaselineG.SetPoint(i,mass,median) 
  oneSigBaselineG.SetPointError(i,0.0,0.0,low1sig,high1sig) 
  twoSigBaselineG.SetPointError(i,0.0,0.0,low2sig,high2sig) 

for i,entry in zip(range(len(limits20Bern4)),limits20Bern4):
  if (float(entry[0])==155):
    continue
  obs20Bern4G.SetPoint(i,float(entry[0]),float(entry[1]))
  exp20Bern4G.SetPoint(i,float(entry[0]),float(entry[4]))

for i,entry in zip(range(len(limits20Bern2)),limits20Bern2):
  if (float(entry[0])==155):
    continue
  obs20Bern2G.SetPoint(i,float(entry[0]),float(entry[1]))
  exp20Bern2G.SetPoint(i,float(entry[0]),float(entry[4]))

for i,entry in zip(range(len(limits20Baseline)),limits20Baseline):
  if (float(entry[0])==155):
    continue
  obs20BaselineG.SetPoint(i,float(entry[0]),float(entry[1]))
  exp20BaselineG.SetPoint(i,float(entry[0]),float(entry[4]))

###############################################################33

obsBaselineG.SetLineColor(root.kBlack)
obsBaselineG.SetMarkerColor(root.kBlack)
expBaselineG.SetLineColor(root.kBlack)
expBaselineG.SetMarkerColor(root.kBlack)

obs20Bern4G.SetLineColor(root.kRed)
obs20Bern4G.SetMarkerColor(root.kRed)
exp20Bern4G.SetLineColor(root.kRed)
exp20Bern4G.SetMarkerColor(root.kRed)

obs20Bern2G.SetLineColor(root.kMagenta)
obs20Bern2G.SetMarkerColor(root.kMagenta)
exp20Bern2G.SetLineColor(root.kMagenta)
exp20Bern2G.SetMarkerColor(root.kMagenta)

obs20BaselineG.SetLineColor(root.kBlue)
obs20BaselineG.SetMarkerColor(root.kBlue)
exp20BaselineG.SetLineColor(root.kBlue)
exp20BaselineG.SetMarkerColor(root.kBlue)

#expBaselineG.SetLineStyle(2)
#exp20Bern4G.SetLineStyle(2)
#exp20BaselineG.SetLineStyle(2)

###############################################################33

canvas = root.TCanvas("canvas")

tlatex = root.TLatex()
tlatex.SetNDC()
tlatex.SetTextFont(root.gStyle.GetLabelFont())
tlatex.SetTextSize(0.04)

axisHist = root.TH2F("axisHist","",1,110,160,1,0,150)
setHistTitles(axisHist,"M(#mu#mu) [GeV/c^{2}]","95% CL Limit on #sigma/#sigma_{SM} (H#rightarrow#mu#mu)")
axisHist.Draw()

#twoSigBaselineG.Draw("3")
#oneSigBaselineG.Draw("3")
exp20Bern4G.Draw("L")
exp20Bern2G.Draw("L")
exp20BaselineG.Draw("L")
expBaselineG.Draw("L")
#obs20Bern4G.Draw("LP")
#obs20BaselineG.Draw("LP")
#obsBaselineG.Draw("LP")

#leg = root.TLegend(0.172,0.905,0.73,0.585)
leg = root.TLegend(0.172,0.905,0.73,0.635)
leg.SetLineColor(0)
leg.SetFillColor(0)
#leg.AddEntry(obsBaselineG,"Observed PAS Result","lp")
#leg.AddEntry(expBaselineG,"Expected PAS Result","l")
#leg.AddEntry(obs20Bern4G,"Observed 20 GeV Window 4-Bern","lp")
#leg.AddEntry(exp20Bern4G,"Expected 20 GeV Window 4-Bern","l")
#leg.AddEntry(obs20BaselineG,"Observed 20 GeV Window Baseline","lp")
#leg.AddEntry(exp20BaselineG,"Expected 20 GeV Window Baseline","l")
leg.AddEntry(expBaselineG,"Baseline (PAS)","l")
leg.AddEntry(exp20BaselineG,"Baseline 20 GeV Window","l")
leg.AddEntry(exp20Bern2G,"2-Bernstein 20 GeV Window","l")
leg.AddEntry(exp20Bern4G,"4-Bernstein 20 GeV Window","l")
leg.Draw()

tlatex.SetTextAlign(12)
tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
tlatex.SetTextAlign(32)
tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,"0,1-Jet Tight BB #sqrt{s}=8 TeV")

saveAs(canvas,"output/limitsCompareWidths")

#raw_input("Enter to Continue...")

