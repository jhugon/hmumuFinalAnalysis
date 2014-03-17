#! /usr/bin/env python

import singleHelpers
import ROOT as root
root.gROOT.SetBatch(True)
import scipy
import scipy.stats

from helpers import *

canvas = root.TCanvas("canvas")

ymax = 2.0

ndf = range(1,101)
medianChi2 = [scipy.stats.chi2.isf(0.5,i) for i in ndf]
low1sigChi2 = [scipy.stats.chi2.isf(scipy.stats.norm.sf(1.),i) for i in ndf]
high1sigChi2 = [scipy.stats.chi2.isf(scipy.stats.norm.sf(-1.),i) for i in ndf]
low1sigChi2 = [abs(i-j) for i,j in zip(low1sigChi2,medianChi2)]
high1sigChi2 = [abs(i-j) for i,j in zip(high1sigChi2,medianChi2)]

medianNormChi2 = [i/j for i,j in zip(medianChi2,ndf)]
low1sigNormChi2 = [scipy.stats.chi2.isf(scipy.stats.norm.sf(-1.),i)/i for i in ndf]
high1sigNormChi2 = [min(scipy.stats.chi2.isf(scipy.stats.norm.sf(1.),i)/i,ymax) for i in ndf]
low1sigNormChi2 = [abs(i-j) for i,j in zip(low1sigNormChi2,medianNormChi2)]
high1sigNormChi2 = [abs(i-j) for i,j in zip(high1sigNormChi2,medianNormChi2)]

low2sigNormChi2 = [scipy.stats.chi2.isf(scipy.stats.norm.sf(-2.),i)/i for i in ndf]
high2sigNormChi2 = [min(scipy.stats.chi2.isf(scipy.stats.norm.sf(2.),i)/i,ymax) for i in ndf]
low2sigNormChi2 = [abs(i-j) for i,j in zip(low2sigNormChi2,medianNormChi2)]
high2sigNormChi2 = [abs(i-j) for i,j in zip(high2sigNormChi2,medianNormChi2)]

chi2Vndf = root.TGraphAsymmErrors()
chi2Vndf.SetLineColor(root.kBlue)
chi2Vndf.SetFillColor(root.kCyan)

normChi2Vndf = root.TGraph()
normChi2Vndf1Sig = root.TGraphAsymmErrors()
normChi2Vndf2Sig = root.TGraphAsymmErrors()

normChi2Vndf.SetLineColor(root.kBlack)
normChi2Vndf1Sig.SetLineColor(root.kGreen)
normChi2Vndf2Sig.SetLineColor(root.kYellow)
normChi2Vndf1Sig.SetFillColor(root.kGreen)
normChi2Vndf2Sig.SetFillColor(root.kYellow)

for i in range(len(ndf)):
  chi2Vndf.SetPoint(i,ndf[i],medianChi2[i])
  chi2Vndf.SetPointError(i,0,0,low1sigChi2[i],high1sigChi2[i])
  normChi2Vndf.SetPoint(i,ndf[i],medianNormChi2[i])
  normChi2Vndf1Sig.SetPoint(i,ndf[i],medianNormChi2[i])
  normChi2Vndf1Sig.SetPointError(i,0,0,low1sigNormChi2[i],high1sigNormChi2[i])
  normChi2Vndf2Sig.SetPoint(i,ndf[i],medianNormChi2[i])
  normChi2Vndf2Sig.SetPointError(i,0,0,low2sigNormChi2[i],high2sigNormChi2[i])

axisHistChi2 = root.TH2F("","",1,0,100,1,0,100)
setHistTitles(axisHistChi2,"NDF","50-Percentile #chi^{2} Value")
axisHistNormChi2 = root.TH2F("","",1,0,100,1,0,ymax)
setHistTitles(axisHistNormChi2,"NDF","#chi^{2}/NDF Value")

axisHistChi2.Draw()
chi2Vndf.Draw("4")
chi2Vndf.Draw("L")
saveAs(canvas,"output/chi2Prob")

canvas.Clear()

axisHistNormChi2.Draw()
normChi2Vndf2Sig.Draw("4")
normChi2Vndf1Sig.Draw("4")
normChi2Vndf.Draw("L")
leg = root.TLegend(0.58,0.70,0.9,0.9)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.AddEntry(normChi2Vndf,"#chi^{2} Dist. Median","l")
leg.AddEntry(normChi2Vndf1Sig,"#chi^{2} Dist. 1#sigma Interval","f")
leg.AddEntry(normChi2Vndf2Sig,"#chi^{2} Dist. 2#sigma Interval","f")
leg.Draw()
canvas.RedrawAxis()
saveAs(canvas,"output/normChi2Prob")
