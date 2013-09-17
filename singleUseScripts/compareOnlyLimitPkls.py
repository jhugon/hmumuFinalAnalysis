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

limits160F = open("limitsCombSplitAll7P8TeV_bugEff.pkl")
limits170F = open("limitsCombSplitAll7P8TeV_fixEff.pkl")

limits160 = cPickle.load(limits160F)
limits170 = cPickle.load(limits170F)

obs160G = root.TGraph()
obs170G = root.TGraph()
exp160G = root.TGraph()
exp170G = root.TGraph()

oneSig160G = root.TGraphAsymmErrors()
oneSig160G.SetFillColor(root.kGreen)
oneSig160G.SetLineStyle(0)
twoSig160G = root.TGraphAsymmErrors()
twoSig160G.SetFillColor(root.kYellow)
twoSig160G.SetLineStyle(0)

for i,entry in zip(range(len(limits160)),limits160):

  mass = float(entry[0])
  obs160G.SetPoint(i,mass,float(entry[1]))

  median = float(entry[4])
  low2sig = median - float(entry[2])
  low1sig = median - float(entry[3])
  high1sig = float(entry[5]) - median
  high2sig = float(entry[6]) - median

  exp160G.SetPoint(i,mass,median)
  oneSig160G.SetPoint(i,mass,median) 
  twoSig160G.SetPoint(i,mass,median) 
  oneSig160G.SetPointError(i,0.0,0.0,low1sig,high1sig) 
  twoSig160G.SetPointError(i,0.0,0.0,low2sig,high2sig) 

  #print 
for i,entry in zip(range(len(limits170)),limits170):
  obs170G.SetPoint(i,float(entry[0]),float(entry[1]))
  exp170G.SetPoint(i,float(entry[0]),float(entry[4]))

###############################################################33

obs160G.SetLineColor(root.kBlue)
obs160G.SetMarkerColor(root.kBlue)
exp160G.SetLineColor(root.kBlue)
exp160G.SetMarkerColor(root.kBlue)

obs170G.SetLineColor(root.kRed)
obs170G.SetMarkerColor(root.kRed)
exp170G.SetLineColor(root.kRed)
exp170G.SetMarkerColor(root.kRed)

exp160G.SetLineStyle(2)
exp170G.SetLineStyle(2)

###############################################################33

canvas = root.TCanvas("canvas")

tlatex = root.TLatex()
tlatex.SetNDC()
tlatex.SetTextFont(root.gStyle.GetLabelFont())
tlatex.SetTextSize(0.04)

axisHist = root.TH2F("axisHist","",1,110,160,1,0,50)
setHistTitles(axisHist,"M(#mu#mu) [GeV/c^{2}]","95% CL Limit on #sigma/#sigma_{SM} (H#rightarrow#mu#mu)")
axisHist.Draw()

twoSig160G.Draw("3")
oneSig160G.Draw("3")
exp170G.Draw("L")
exp160G.Draw("L")
obs170G.Draw("LP")
obs160G.Draw("LP")

leg = root.TLegend(0.172,0.905,0.73,0.585)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.AddEntry(obs160G,"Observed M(#mu#mu), efficiency #in [110,160] GeV/c^{2}","lp")
leg.AddEntry(exp160G,"Expected M(#mu#mu), efficiency #in [110,160] GeV/c^{2}","l")
leg.AddEntry(obs170G,"Observed M(#mu#mu), efficiency #in [110,170] GeV/c^{2}","lp")
leg.AddEntry(exp170G,"Expected M(#mu#mu), efficiency #in [110,170] GeV/c^{2}","l")
leg.Draw()

tlatex.SetTextAlign(12)
tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
tlatex.SetTextAlign(32)
tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,"H#rightarrow#mu#mu Combination 2011+2012")

saveAs(canvas,"output/limitsCompare")

