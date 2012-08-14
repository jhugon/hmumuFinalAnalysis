#!/usr/bin/python

import math
import ROOT as root
from helpers import *
import matplotlib.pyplot as mpl
import numpy

from xsec import *

directory = "input/open/"

analysisList = ["","VBFSelected","VBFTightSelected","ZPt30Selected","ZPt75Selected"]

signalName = "ggHmumu125"
backgroundName = "DYJetsToLL"

scaleSignal = 1.0

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

class DataFromHists:
  def __init__(self,directory,signalName,backgroundName,analysisList,massRange=[123,127]):
    self.sigFile = root.TFile(directory+signalName+".root")
    self.bakFile = root.TFile(directory+backgroundName+".root")
    analyses = analysisList
    self.analyses = analyses
    self.sigHists = []
    self.bakHists = []
    for a in analyses:
        self.sigHists.append(self.sigFile.Get("mDiMu"+a))
        self.bakHists.append(self.bakFile.Get("mDiMu"+a))

    effMap = {}
    xsecMap = {}
    for sigHist,bakHist,anaName in zip(self.sigHists,self.bakHists,analyses):
      axis = sigHist.GetXaxis()
      lowBin = axis.FindBin(massRange[0])
      assert(axis.GetBinLowEdge(lowBin)==massRange[0])
      highBin = axis.FindBin(massRange[1])
      highBin -= 1
      assert(axis.GetBinUpEdge(highBin)==massRange[1])
      countsSig = sigHist.Integral(lowBin,highBin)
      countsBak = bakHist.Integral(lowBin,highBin)
      
      effSig = countsSig/nEventsMap[signalName]
      effBak = countsBak/nEventsMap[backgroundName]
      xsecSig = effSig*xsec[signalName]
      xsecBak = effBak*xsec[backgroundName]
      effMap[anaName] = {"signal":effSig,"background":effBak}
      xsecMap[anaName] = {"signal":xsecSig,"background":xsecBak}
    self.effMap = effMap
    self.xsecMap = xsecMap

  def getSigEff(self,anaName):
    return self.effMap[anaName]["signal"]
  def getBakXSec(self,anaName):
    return self.xsecMap[anaName]["background"]

  def calc(self,lumi,anaName,scaleSignal=1.0):
    lumi = lumi*1000.0
    s = self.xsecMap[anaName]["signal"]*lumi*scaleSignal
    b = self.xsecMap[anaName]["background"]*lumi
    #fom = s/(1.5+math.sqrt(b)+self.backUnc*b)
    #sig = s/sqrt(b+s)
    sig = s/sqrt(b)
    return s, b, sig

###### Text Part

textfile = open("nEvents.txt","w")
texfile = open("nEvents.tex","w")

#Different Mass Ranges
textfile.write("\nTesting Different Search Windows, Inclusive Analysis\n")
textfile.write("{0:<13} {1:<13} {2:<13}\n".format("125+/-X [GeV]","Sig Eff","Bak Xsec [pb]"))
for i in range(10):
  i = 0.5+0.5*i
  data = DataFromHists(directory,signalName,backgroundName,[""],massRange=[125.-i,125.+i])
  eff = data.getSigEff("")
  back = data.getBakXSec("")
  textfile.write("{0:<13.2f} {1:<13.3e} {2:<13.3e}\n".format(i,eff,back))

textfile.write("\n#######################################\n")

data = DataFromHists(directory,signalName,backgroundName,analysisList)

for ana in analysisList:
  textfile.write("\n{0} {1} {2}\n".format(ana,signalName,backgroundName))
  textfile.write("{0:<13} {1:<13}\n".format("Sig Eff","Bak Xsec [pb]"))
  eff = data.getSigEff(ana)
  back = data.getBakXSec(ana)
  textfile.write("{0:<13.3e} {1:<13.3e}\n".format(eff,back))

textfile.write("\n#######################################\n")
  
for lumi in lumis:
  textfile.write("\nIntegrated Luminosity (fb^-1): {0:.2f}\n".format(lumi))
  texfile.write("\\begin{tabular{|l|l|l|l|l|l|}} \\hline \n")
  texfile.write("\\%multicolumn{0}{{Events For {1:.2f}fb$^{2}$}} \\\\ \\hline\n".format("{6}",lumi,"{-1}"))
  texfile.write("Selection & Signal & $S$ & $B$ & $S/B$ & $S/\\sqrt{B}$ \\\\ \\hline \n")
  if scaleSignal != 1.0:
    textfile.write("\nHiggs Signal Scaled by factor of: {0:.2f}\n".format(scaleSignal))
  for ana in analysisList:
      s, b, sosqrtb = data.calc(lumi,ana,scaleSignal=scaleSignal)
      textfile.write("{0:<17} {1:<10} s: {2:<6.2f} b: {3:<10.2f} s/b: {4:<6.2e} s/sqrt(b): {5:<5.2f}\n".format(ana,signalName,s,b,s/b,sosqrtb))
      texfile.write("{0:<17} & {1:<10} & {2:<6.2f} & {3:<10.2f} & {4:<6.2e} & {5:<5.2f} \\\\ \\hline \n".format(ana,signalName,s,b,s/b,sosqrtb))

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
for ana in analysisList:
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
  for ana in analysisList:
      s, b, sosqrtb = data.calc(lumi,ana)
      significancePlots[ana].SetPoint(iPoint,lumi,sosqrtb)
  iPoint += 1

iPoint = 0
for lumi in lumis:
  for ana in analysisList:
      found3Sig = False
      found5Sig = False
      #rangeToUse = numpy.logspace(1.0,100.0,num=1000)
      rangeToUse = [0.1,0.5,1.0,2.0,5.0,8.,10.,12.,15.,20.,30.,40.,50.,75.,150.,200,300.,400.,500.,700.,1000.,1500.,2000.,5000.]
      for mult in rangeToUse:
        s, b, sosqrtb = data.calc(lumi,ana,scaleSignal=mult)
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
for ana in analysisList:
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
for ana in analysisList:
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
for ana in analysisList:
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


"""
sobList = []
for ana in analysisList:
      s, b, sosqrtb = data.calc(lumi,ana)
      sobList.append(s/b)
      break

fig = mpl.figure()
ax1 = fig.add_subplot(111)
ax1bounds = ax1.get_position().bounds
ax1.set_position([0.2,0.1,0.7,0.85]) #uncomment if you need more space for names
#ax1.set_position([0.25,0.1,0.7,0.85]) #uncomment if you need more space for names
pos = numpy.arange(len(analysisList))
ax1.grid(axis="x")
ax1.set_yticks(pos+0.25)
ax1.set_yticklabels(tuple(analysisList))
ax1.set_xlabel("S/B")
bars = ax1.barh(pos,sobList, 0.5,log=True)
fig.savefig("sob.png")
fig.savefig("sob.pdf")
"""

