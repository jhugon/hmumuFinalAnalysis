#!/usr/bin/python

import math
import ROOT as root
from helpers import *

xsec = {}
xsec["ggHmumu"] = 4.236e-3
xsec["vbfHmumu"] = 3.338e-4
xsec["WHmumu"] = 1.512e-4
xsec["ZHmumu"] = 8.556e-5

eff = {}
eff["inc"] = {
  "ggHmumu": 0.635
}
eff["vbf"] = {
  "vbfHmumu": 0.115
}
eff["vbfTight"] = {
  "vbfHmumu": 0.0664
}
eff["pt50"] = {
  "ggHmumu": 0.185
}
eff["pt75"] = {
  "ggHmumu": 0.106
}
eff["incBZMuCuts"] = {
  "ggHmumu": 0.397
}
eff["pt50BZMuCuts"] = {
  "ggHmumu": 0.0789
}

back = {}
back["inc"] = 0.8
back["vbf"] = 8e-4
back["vbfTight"] = 1.6e-4
back["pt50"] = 0.1
back["pt75"] = 4.8e-2
back["incBZMuCuts"] = 0.36
back["pt50BZMuCuts"] = 0.032

lumis = [
5.0,
10.0,
20.0,
30.0,
50.#,
#75.,
#100.,
#200.,
#500.,
#1000.
]

colors = [
root.kBlue+1,
root.kRed+1,
root.kGreen+1,
root.kOrange,
1,
root.kCyan,
root.kPink
]

###### Text Part

textfile = open("nEvents.txt","w")

scaleSignal = 1.0

class Calculator:
  def __init__(self,xsec,eff,back):
    self.xsec = xsec
    self.eff = eff
    self.back = back

  def calc(self,lumi,analysis,signal,scaleSignal=1.0):
    lumi = lumi*1000.0
    s = self.xsec[signal]*self.eff[analysis][signal]*lumi*scaleSignal
    b = self.back[analysis]*lumi
    return s, b, s/math.sqrt(b)
  
analysis = Calculator(xsec,eff,back)

sortedAnalysis = eff.keys()
sortedAnalysis.sort()
for lumi in lumis:
  textfile.write("\nIntegrated Luminosity (fb^-1): {0:.2f}\n".format(lumi))
  if scaleSignal != 1.0:
    textfile.write("\nHiggs Signal Scaled by factor of: {0:.2f}\n".format(scaleSignal))
  for ana in sortedAnalysis:
    for signal in eff[ana]:
      s, b, sosqrtb = analysis.calc(lumi,ana,signal,scaleSignal=scaleSignal)
      textfile.write("{0:<15} {1:<10} s: {2:<6.2f} b: {3:<10.2f} s/b: {4:<6.2e} s/sqrt(b): {5:<5.2f}\n".format(ana,signal,s,b,s/b,sosqrtb))

textfile.close()

########## End Text Part

root.gROOT.SetBatch(True)
setStyle()

canvas = root.TCanvas("canvas")
#canvas.SetLogx(1)
#canvas.SetLogy(1)

leg = root.TLegend(0.75,0.75,0.9,0.9)
leg.SetFillColor(0)
leg.SetLineColor(0)

significancePlots = {}
sobPlots = {}
iColor = 0
for ana in sortedAnalysis:
  tmp = root.TGraph()
  tmp.SetLineColor(colors[iColor % len(colors)])
  significancePlots[ana] = tmp

  leg.AddEntry(tmp,ana,"l")

  tmp2 = root.TGraph()
  tmp2.SetLineColor(colors[iColor % len(colors)])
  sobPlots[ana] = tmp2
  iColor += 1

iPoint = 0
for lumi in lumis:
  for ana in sortedAnalysis:
    for signal in eff[ana]:
      s, b, sosqrtb = analysis.calc(lumi,ana,signal)
      significancePlots[ana].SetPoint(iPoint,lumi,sosqrtb)
      sobPlots[ana].SetPoint(iPoint,lumi,s/b)
      print s/b
  iPoint += 1

firstPlot = True
for ana in sortedAnalysis:
  if firstPlot:
    firstPlot = False
    significancePlots[ana].Draw("al")
    significancePlots[ana].GetXaxis().SetTitle("Integrated Luminosity [fb^{-1}]")
    significancePlots[ana].GetYaxis().SetTitle("S/#sqrt{B}")
    significancePlots[ana].Draw("al")

  else:
    significancePlots[ana].Draw("l")
leg.Draw()
saveAs(canvas,"significance")

canvas.SetLogy(1)
firstPlot = True
for ana in sortedAnalysis:
  if firstPlot:
    firstPlot = False
    sobPlots[ana].Draw("al")
    sobPlots[ana].GetXaxis().SetTitle("Integrated Luminosity [fb^{-1}]")
    sobPlots[ana].GetYaxis().SetTitle("S/B")
    sobPlots[ana].GetYaxis().SetRangeUser(0.001,0.2)
    sobPlots[ana].Draw("al")

  else:
    sobPlots[ana].Draw("l")
leg.Draw()
saveAs(canvas,"sob")
