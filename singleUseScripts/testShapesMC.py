#! /usr/bin/env python

from ROOT import gSystem

import sys
import os
import re
import math
from ROOT import *
gSystem.Load('libRooFit')
import ROOT as root

import singleHelpers
from helpers import *
from xsec import *
import makeCards
import fitOrderChooser

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT

#####################################################################

maxMass = 160.
minMass = 110.

category = "Jets01PassPtG10BB"
categoryTitle = "0,1-Jet Tight BB"
energyStr = "8TeV"

order = None

pdfFunc = makeCards.makePDFBakExpMOverSq
#pdfFunc = makeCards.makePDFBakOld
#pdfFunc = fitOrderChooser.makePDFBakBernstein
#pdfFunc = fitOrderChooser.makePDFBakSumExp

directory = "/data/uftrig01b/jhugon/hmumu/analysisV00-01-10/forGPReRecoMuScleFit/"
directory = "/afs/cern.ch/work/j/jhugon/public/hmumuNtuplesLevel2/unzipped/"

useWeight2Error = True

#####################################################################

fDY = root.TFile(directory+"DYJetsToLL_"+energyStr+".root")
fTT = root.TFile(directory+"ttbar_"+energyStr+".root")

tDY = fDY.Get("outtree"+category)
tTT = fTT.Get("outtree"+category)

#####################################################################

controlRegionVeryLow=[60,110]
controlRegionLow=[110,120]
controlRegionHigh=[130,160]
dimuonMass = root.RooRealVar("dimuonMass","M(#mu#mu) [GeV/c^{2}]",minMass,maxMass)
dimuonMass.setRange("exprange",120.,controlRegionHigh[1])
dimuonMass.setRange("whole",controlRegionLow[0],controlRegionHigh[1])
dimuonMass.setRange("low",controlRegionLow[0],controlRegionLow[1])
dimuonMass.setRange("high",controlRegionHigh[0],controlRegionHigh[1])
dimuonMass.setRange("signal",controlRegionLow[1],controlRegionHigh[0])
dimuonMass.setRange("signalfit",110,140)
dimuonMass.setBins(50)

weightDY = lumiDict[energyStr]*1000.*xsec["DYJetsToLL_"+energyStr]/nEventsMap["DYJetsToLL_"+energyStr]
weightTT = lumiDict[energyStr]*1000.*xsec["ttbar_"+energyStr]/nEventsMap["ttbar_"+energyStr]
print "weightDY: ",weightDY
print "weightTT: ",weightTT

weightVar = root.RooRealVar("weight","Event Weight",1.)

dataNoWeightDY = root.RooDataSet("dataNoWeightDY","No Weight DY Data",tDY,root.RooArgSet(dimuonMass))
dataNoWeightTT = root.RooDataSet("dataNoWeightTT","No Weight TT Data",tTT,root.RooArgSet(dimuonMass))

dataNoWeightDY.Print()
dataNoWeightTT.Print()

weightVar.setVal(weightDY)
dataDY = root.RooDataSet("dataDY","DY Data",root.RooArgSet(dimuonMass,weightVar),root.RooFit.Import(dataNoWeightDY),root.RooFit.WeightVar(weightVar))
weightVar.setVal(weightTT)
dataTT = root.RooDataSet("dataTT","TT Data",root.RooArgSet(dimuonMass,weightVar),root.RooFit.Import(dataNoWeightTT),root.RooFit.WeightVar(weightVar))
weightVar.setVal(weightDY)
dataBoth = root.RooDataSet("dataBoth","DY+TT Data",dataDY,root.RooArgSet(dimuonMass,weightVar),"","weight")
weightVar.setVal(weightTT)
dataBoth.append(dataTT)

dataDY.Print()
print dataDY.sumEntries()
dataTT.Print()
print dataTT.sumEntries()
dataBoth.Print()
print dataBoth.sumEntries()

#####################################################################

w = root.RooWorkspace("w")
wImport = getattr(w,"import")

dimuonMassZ = None
dataZ = None

tmpParamList,tmpNormTup,tmpDebug,tmpOrder = pdfFunc(category+energyStr,dataDY,dimuonMass,minMass,maxMass,wImport,dimuonMassZ,dataZ,order=order)
pdf = w.pdf("bak")

#####################################################################

canvas = root.TCanvas("canvas")

frDY = pdf.fitTo(dataDY,root.RooFit.Save(),PRINTLEVEL,root.RooFit.SumW2Error(useWeight2Error))
frDY.Print()
rmpDY = RooModelPlotter(dimuonMass,pdf,dataDY,frDY,
                      categoryTitle,energyStr,lumiDict[energyStr],
                      caption2="Drell-Yan MC",
                      canvas=canvas
                      )
rmpDY.draw("mcFit_"+category+"DY")
canvas.Clear()

frTT = pdf.fitTo(dataTT,root.RooFit.Save(),PRINTLEVEL,root.RooFit.SumW2Error(useWeight2Error))
frTT.Print()
rmpTT = RooModelPlotter(dimuonMass,pdf,dataTT,frTT,
                      categoryTitle,energyStr,lumiDict[energyStr],
                      caption2=r"t\bar{t} MC",
                      canvas=canvas
                      )
rmpTT.draw("mcFit_"+category+"TT")

frBoth = pdf.fitTo(dataBoth,root.RooFit.Save(),PRINTLEVEL,root.RooFit.SumW2Error(useWeight2Error))
frBoth.Print()
rmpBoth = RooModelPlotter(dimuonMass,pdf,dataBoth,frBoth,
                      categoryTitle,energyStr,lumiDict[energyStr],
                      caption2=r"DY+t\bar{t} MC",
                      showStackDatasets=[dataTT,dataDY],
                      showStackDatasetTitles=[legendEntries['ttbar'],legendEntries["DYJetsToLL"]],
                      showStackDatasetColors=[colors['ttbar'],colors["DYJetsToLL"]],
                      canvas=canvas
                      )
rmpBoth.draw("mcFit_"+category+"Both")
