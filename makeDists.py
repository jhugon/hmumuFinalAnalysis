#!/usr/bin/python

import math
import ROOT as root
from helpers import *
import matplotlib.pyplot as mpl
import numpy

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
eff["vbfLoose"] = {
  "vbfHmumu": 0.173
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
back["vbfLoose"] = 1.4e-3
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
root.kPink,
root.kMagenta,
root.kSpring-9
]

BACKUNC=0.15

###### Text Part

textfile = open("nEvents.txt","w")
texfile = open("nEvents.tex","w")

scaleSignal = 1.0

class Calculator:
  def __init__(self,xsec,eff,back,backUnc):
    self.xsec = xsec
    self.eff = eff
    self.back = back
    self.backUnc = backUnc

  def calc(self,lumi,analysis,signal,scaleSignal=1.0):
    lumi = lumi*1000.0
    s = self.xsec[signal]*self.eff[analysis][signal]*lumi*scaleSignal
    b = self.back[analysis]*lumi
    #fom = s/(1.5+math.sqrt(b)+self.backUnc*b)
    #sig = s/sqrt(b+s)
    sig = s/sqrt(b)
    return s, b, sig
  
analysis = Calculator(xsec,eff,back,BACKUNC)

sortedAnalysis = eff.keys()
sortedAnalysis.sort()
for lumi in lumis:
  textfile.write("\nIntegrated Luminosity (fb^-1): {0:.2f}\n".format(lumi))
  texfile.write("\\begin{tabular{|l|l|l|l|l|l|}} \\hline \n")
  texfile.write("\\%multicolumn{0}{{Events For {1:.2f}fb$^{2}$}} \\\\ \\hline\n".format("{6}",lumi,"{-1}"))
  texfile.write("Selection & Signal & $S$ & $B$ & $S/B$ & $S/\\sqrt{B}$ \\\\ \\hline \n")
  if scaleSignal != 1.0:
    textfile.write("\nHiggs Signal Scaled by factor of: {0:.2f}\n".format(scaleSignal))
  for ana in sortedAnalysis:
    for signal in eff[ana]:
      s, b, sosqrtb = analysis.calc(lumi,ana,signal,scaleSignal=scaleSignal)
      textfile.write("{0:<15} {1:<10} s: {2:<6.2f} b: {3:<10.2f} s/b: {4:<6.2e} s/sqrt(b): {5:<5.2f}\n".format(ana,signal,s,b,s/b,sosqrtb))
      texfile.write("{0:<15} & {1:<10} & {2:<6.2f} & {3:<10.2f} & {4:<6.2e} & {5:<5.2f} \\\\ \\hline \n".format(ana,signal,s,b,s/b,sosqrtb))

  texfile.write("\\hline\n\\end{tabular}\n\n")

textfile.close()
texfile.close()

########## End Text Part

root.gROOT.SetBatch(True)
setStyle()

canvas = root.TCanvas("canvas")
#canvas.SetLogx(1)
#canvas.SetLogy(1)

#leg = root.TLegend(0.75,0.75,0.9,0.9)
leg = root.TLegend(0.23,0.6,0.5,0.9)
leg.SetFillColor(0)
leg.SetLineColor(0)

significancePlots = {}
threeSigPlots = {}
fiveSigPlots = {}
iColor = 0
for ana in sortedAnalysis:
  tmp = root.TGraph()
  tmp.SetLineColor(colors[iColor % len(colors)])
  significancePlots[ana] = tmp

  tmp = root.TGraph()
  tmp.SetLineColor(colors[iColor % len(colors)])
  threeSigPlots[ana] = tmp

  tmp = root.TGraph()
  tmp.SetLineColor(colors[iColor % len(colors)])
  fiveSigPlots[ana] = tmp


  leg.AddEntry(tmp,ana,"l")

  iColor += 1

iPoint = 0
for lumi in lumis:
  for ana in sortedAnalysis:
    for signal in eff[ana]:
      s, b, sosqrtb = analysis.calc(lumi,ana,signal)
      significancePlots[ana].SetPoint(iPoint,lumi,sosqrtb)
  iPoint += 1

iPoint = 0
for lumi in lumis:
  for ana in sortedAnalysis:
    for signal in eff[ana]:
      found3Sig = False
      found5Sig = False
      #rangeToUse = numpy.logspace(1.0,100.0,num=1000)
      rangeToUse = [0.1,0.5,1.0,2.0,5.0,8.,10.,12.,15.,20.,30.,40.,50.,75.,150.,200,300.,400.,500.,700.,1000.,1500.,2000.,5000.]
      for mult in rangeToUse:
        s, b, sosqrtb = analysis.calc(lumi,ana,signal,scaleSignal=mult)
        if sosqrtb>3.0 and not found3Sig:
          threeSigPlots[ana].SetPoint(iPoint,lumi,mult)
          found3Sig = True
        if sosqrtb>5.0:
          fiveSigPlots[ana].SetPoint(iPoint,lumi,mult)
          found5Sig = True
      if not found3Sig:
          threeSigPlots[ana].SetPoint(iPoint,lumi,1000.0)
      if not found5Sig:
          fiveSigPlots[ana].SetPoint(iPoint,lumi,10000.0)
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
    threeSigPlots[ana].Draw("al")
    threeSigPlots[ana].GetXaxis().SetTitle("Integrated Luminosity [fb^{-1}]")
    threeSigPlots[ana].GetYaxis().SetTitle("S/#sqrt{B} #geq 3 for SM #times X")
    threeSigPlots[ana].GetYaxis().SetRangeUser(1,100)
    threeSigPlots[ana].Draw("al")

  else:
    threeSigPlots[ana].Draw("l")
leg.Draw()
saveAs(canvas,"threeSig")

canvas.SetLogy(1)
firstPlot = True
for ana in sortedAnalysis:
  if firstPlot:
    firstPlot = False
    fiveSigPlots[ana].Draw("al")
    fiveSigPlots[ana].GetXaxis().SetTitle("Integrated Luminosity [fb^{-1}]")
    fiveSigPlots[ana].GetYaxis().SetTitle("S/#sqrt{B} #geq 5 for SM #times X")
    fiveSigPlots[ana].GetYaxis().SetRangeUser(100,10000)
    fiveSigPlots[ana].Draw("al")

  else:
    fiveSigPlots[ana].Draw("l")
leg.Draw()
saveAs(canvas,"fiveSig")


sobList = []
for ana in sortedAnalysis:
  for signal in eff[ana]:
      s, b, sosqrtb = analysis.calc(lumi,ana,signal)
      sobList.append(s/b)
      break

fig = mpl.figure()
ax1 = fig.add_subplot(111)
ax1bounds = ax1.get_position().bounds
ax1.set_position([0.2,0.1,0.7,0.85]) #uncomment if you need more space for names
#ax1.set_position([0.25,0.1,0.7,0.85]) #uncomment if you need more space for names
pos = numpy.arange(len(sortedAnalysis))
ax1.grid(axis="x")
ax1.set_yticks(pos+0.25)
ax1.set_yticklabels(tuple(sortedAnalysis))
ax1.set_xlabel("S/B")
bars = ax1.barh(pos,sobList, 0.5,log=True)
fig.savefig("sob.png")
fig.savefig("sob.pdf")
