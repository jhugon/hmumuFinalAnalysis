#!/usr/bin/env python

from helpers import *
import math
import os.path
import glob
import numpy
import matplotlib.pyplot as mpl

root.gErrorIgnoreLevel = root.kWarning
root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

infilename = "input/ggHmumu125_8TeV.root"
infilename = "input/smearing/ggHmumu125_8TeV.root"
infilename = "input/smearing/vbfHmumu125_8TeV.root"

categories = ["BB","BO","BE","OO","OE","EE"]
keyList = categories + ["NotBB","All"]
histMap = {}
countsMap = {}

colors = [root.kRed,root.kBlue,root.kGreen,root.kOrange,root.kPink,root.kCyan,root.kBlue+1,root.kBlack]

infile = root.TFile(infilename)
histBase = "mDiMu"

histMap["All"] = infile.Get(histBase)
histMap["NotBB"] = infile.Get("NotBB"+"/"+histBase)
for c in categories:
  histMap[c] = infile.Get(c+"/"+histBase)

for k in histMap:
  countsMap[k] = histMap[k].Integral()

print("{0:<8}{1:}".format("","Fraction Of Signal"))
countsAll = countsMap["All"]
for k in keyList:
  #countsMap[k] = (countsMap[k]) / countsAll
  countsMap[k] /= countsAll
  print("{0:<8}{1:>6.1%}".format(k,countsMap[k]))
print

for c in histMap:
  histMap[c].Scale(1.0/histMap[c].Integral())

print("{0:<8}{1:<10}{2:<6}".format("","Mean","RMS"))
for c in keyList:
  print("{0:<8}{1:<10.2f}{2:<6.2f}".format(c,histMap[c].GetMean(),histMap[c].GetRMS()))
print

resMap={}
print("{0:<8}{1:<10}{2:^24}{3:^24}".format("","Median","1 Sigma Quantile Width","2 Sigma Quantile Width"))
for c in keyList:
  quantiles = getMedianAndQuantileInterval(histMap[c],0.159)
  err = (quantiles[2]-quantiles[0])/2.
  resMap[c] = err
  quantiles2 = getMedianAndQuantileInterval(histMap[c],0.023)
  err2 = (quantiles2[2]-quantiles2[0])/2.
  print("{0:<8}{1:<10.2f}{2:^24.2f}{3:^24.2f}".format(c,quantiles[1],err,err2))
print

for k,c in zip(keyList,colors):
  histMap[k].SetFillStyle(0)
  histMap[k].UseCurrentStyle()
  histMap[k].SetLineColor(c)

drawn = False
canvas = root.TCanvas()
keyList.reverse()
for k in keyList:
  if drawn:
    histMap[k].Draw("hist same")
  else:
    histMap[k].SetTitle("")
    histMap[k].GetXaxis().SetRangeUser(100,150)
    histMap[k].GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    histMap[k].GetYaxis().SetTitle("Normalized Events")
    histMap[k].Draw("hist")
    drawn = True

urLegendPos = [0.70,0.67,0.9,0.9]
leg = root.TLegend(*urLegendPos)
leg.SetFillColor(0)
leg.SetLineColor(0)
for k in keyList:
  leg.AddEntry(histMap[k],k,"l")

leg.Draw()
saveAs(canvas,"resolutionCategories")

##################################################################

fig = mpl.figure()
ax = fig.add_subplot(111)
ax.pie([countsMap[i] for i in categories],labels=tuple(categories),
        shadow=True)
ax.set_title(r"$gg \rightarrow H \rightarrow \mu\mu$ Fractions")
fig.savefig("resFractions.png")

fig.clf()
ax = fig.add_subplot(111)
xPos = numpy.arange(len(categories))+0.25
ax.bar(xPos,[resMap[i] for i in categories],0.5)
ax.set_xticks(xPos+0.25)
ax.set_xticklabels(tuple(categories))
ax.set_title(r"$gg \rightarrow H \rightarrow \mu\mu$ Fractions")
ax.set_ylabel(r"1 $\sigma$ Quantile Width")
fig.savefig("resPlot.png")
