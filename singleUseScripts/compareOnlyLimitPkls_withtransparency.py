#!/usr/bin/env python

from singleHelpers import *
from helpers import *
import ROOT as root
import cPickle
from xsec import *
from scipy.stats import norm

PRELIMINARYSTRING="CMS Preliminary"
root.gStyle.SetMarkerSize(0.9)

###############################################################

#limitsOneF = open("pklfiles/limitsCombSplitAll7P8TeV_voigtExp_newparam.pkl")
limitsOneF = open("pklfiles/limitsCombSplitAll7P8TeV_default.pkl")
limitsTwoF = open("pklfiles/limitsCombSplitAll7P8TeV_voigtexp_zmass_zSig_Constr.pkl")

limitsOne = cPickle.load(limitsOneF)
limitsTwo = cPickle.load(limitsTwoF)

obsOneG = root.TGraph()
obsTwoG = root.TGraph()
expOneG = root.TGraph()
expTwoG = root.TGraph()

oneSigOneG = root.TGraphAsymmErrors()
oneSigOneG.SetFillColor(root.kGreen)
oneSigOneG.SetLineStyle(0)
twoSigOneG = root.TGraphAsymmErrors()
twoSigOneG.SetFillColor(root.kYellow)
twoSigOneG.SetLineStyle(0)

for i,entry in zip(range(len(limitsOne)),limitsOne):

  mass = float(entry[0])
  #print mass
  #if (mass==155):
  #  continue
  obsOneG.SetPoint(i,mass,float(entry[1]))

  median = float(entry[4])
  low2sig = median - float(entry[2])
  low1sig = median - float(entry[3])
  high1sig = float(entry[5]) - median
  high2sig = float(entry[6]) - median

  expOneG.SetPoint(i,mass,median)
  oneSigOneG.SetPoint(i,mass,median) 
  twoSigOneG.SetPoint(i,mass,median) 
  oneSigOneG.SetPointError(i,0.0,0.0,low1sig,high1sig) 
  twoSigOneG.SetPointError(i,0.0,0.0,low2sig,high2sig) 

  #print 
for i,entry in zip(range(len(limitsTwo)),limitsTwo):
  #if (float(entry[0])==155):
  #  continue
  obsTwoG.SetPoint(i,float(entry[0]),float(entry[1]))
  expTwoG.SetPoint(i,float(entry[0]),float(entry[4]))

###############################################################

obsOneG.SetLineColor(root.kRed)
obsOneG.SetMarkerColor(root.kRed)
expOneG.SetLineColor(root.kRed)
expOneG.SetMarkerColor(root.kRed)

obsTwoG.SetLineColor(root.kBlue)
obsTwoG.SetMarkerColor(root.kBlue)
expTwoG.SetLineColor(root.kBlue)
expTwoG.SetMarkerColor(root.kBlue)

expOneG.SetLineStyle(2)
expTwoG.SetLineStyle(2)

###############################################################

canvas = root.TCanvas("canvas")

tlatex = root.TLatex()
tlatex.SetNDC()
tlatex.SetTextFont(root.gStyle.GetLabelFont())
tlatex.SetTextSize(0.04)

axisHist = root.TH2F("axisHist","",1,115,155,1,0,82)
setHistTitles(axisHist,"M_{H} [GeV/c^{2}]","95% CL Limit on #sigma/#sigma_{SM} (H#rightarrow#mu#mu)")
axisHist.Draw()

twoSigOneG.Draw("3")
oneSigOneG.Draw("3")
expOneG.Draw("L")
obsOneG.Draw("LP")
obsTwoG.Draw("LP")

bandHist1 = root.TH1F("bandHist1","",40,115,155)
bandHist1.SetFillStyle(3004)
bandHist1.SetLineWidth(0)
bandHist1.SetLineColor(root.kGray+3)
bandHist1.SetFillColor(root.kGray+3)

bandHist2 = root.TH1F("bandHist2","",40,115,155)
bandHist2.SetFillStyle(3005)
bandHist2.SetLineWidth(0)
bandHist2.SetLineColor(root.kGray+3)
bandHist2.SetFillColor(root.kGray+3)

for i in range(1,6):
  bandHist1.SetBinContent(i,82)
  bandHist2.SetBinContent(i,82)

for i in range(36,41):
  bandHist1.SetBinContent(i,82)
  bandHist2.SetBinContent(i,82)

bandHist1.Draw("same")
bandHist2.Draw("same")
axisHist .Draw("axis same")


leg = root.TLegend(0.33,0.53,0.76,0.73)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.AddEntry(obsOneG,"e^{p_{1} m} / (m-p_{2})^{2}","lp")
leg.AddEntry(obsTwoG,"Voigtian + Exponential","lp")
leg.Draw()

tlatex.SetTextAlign(12)
tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
tlatex.SetTextAlign(32)
tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,"Combination")

saveAs(canvas,"output/limitsComparison_m110to160")

