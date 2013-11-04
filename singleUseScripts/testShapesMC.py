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

TITLEMAP = {
  "Jets01PassPtG10BB": "0,1-Jet Tight BB",
  "Jets01PassPtG10BO": "0,1-Jet Tight BO",
  "Jets01PassPtG10BE": "0,1-Jet Tight BE",
  "Jets01PassPtG10OO": "0,1-Jet Tight OO",
  "Jets01PassPtG10OE": "0,1-Jet Tight OE",
  "Jets01PassPtG10EE": "0,1-Jet Tight EE",
  #"Jets01PassCatAll" : "0,1-Jet Tight Combination",
                        
  "Jets01FailPtG10BB": "0,1-Jet Loose BB",
  "Jets01FailPtG10BO": "0,1-Jet Loose BO",
  "Jets01FailPtG10BE": "0,1-Jet Loose BE",
  "Jets01FailPtG10OO": "0,1-Jet Loose OO",
  "Jets01FailPtG10OE": "0,1-Jet Loose OE",
  "Jets01FailPtG10EE": "0,1-Jet Loose EE",
  #"Jets01FailCatAll" : "0,1-Jet Loose Combination",

  "Jet2CutsVBFPass":"2-Jet VBF Tight",
  "Jet2CutsGFPass":"2-Jet GF Tight",
  "Jet2CutsFailVBFGF":"2-Jet Loose",
}

#####################################################################

maxMass = 160.
minMass = 110.

categoryList = ["Jets01PassPtG10BB","Jets01FailPtG10BB"]
categoryList = TITLEMAP.keys()
categoryList = ["Jet2CutsFailVBFGF"]
energyStr = "8TeV"

order = None

pdfFunc = makeCards.makePDFBakExpMOverSq
#pdfFunc = makeCards.makePDFBakOld
#pdfFunc = fitOrderChooser.makePDFBakBernstein
#pdfFunc = fitOrderChooser.makePDFBakSumExp

directory = "/data/uftrig01b/jhugon/hmumu/analysisV00-01-10/forGPReRecoMuScleFit/"
directory = "/afs/cern.ch/work/j/jhugon/public/hmumuNtuplesLevel2/unzipped/"

useWeight2Error = True
#errorsToDraw = root.RooFit.DataError(root.RooAbsData.Poisson)
errorsToDraw = root.RooFit.DataError(root.RooAbsData.SumW2)


#####################################################################

fDY = root.TFile(directory+"DYJetsToLL_"+energyStr+".root")
fTT = root.TFile(directory+"ttbar_"+energyStr+".root")

for category in categoryList:

  categoryTitle = TITLEMAP[category]

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
  dataTT.Print()
  dataBoth.Print()
  
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
#  frDY.Print()
  rmpDY = RooModelPlotter(dimuonMass,pdf,dataDY,frDY,
                        categoryTitle,energyStr,lumiDict[energyStr],
                        caption2="Drell-Yan MC",
                        errorsToDraw=errorsToDraw,
                        canvas=canvas
                        )
  rmpDY.draw("mcFit_"+category+"_DY")
  canvas.Clear()
  
  frTT = pdf.fitTo(dataTT,root.RooFit.Save(),PRINTLEVEL,root.RooFit.SumW2Error(useWeight2Error))
#  frTT.Print()
  rmpTT = RooModelPlotter(dimuonMass,pdf,dataTT,frTT,
                        categoryTitle,energyStr,lumiDict[energyStr],
                        caption2=r"t\bar{t} MC",
                        errorsToDraw=errorsToDraw,
                        canvas=canvas
                        )
  rmpTT.draw("mcFit_"+category+"_TT")
  
  frBoth = pdf.fitTo(dataBoth,root.RooFit.Save(),PRINTLEVEL,root.RooFit.SumW2Error(useWeight2Error))
#  frBoth.Print()
  rmpBoth = RooModelPlotter(dimuonMass,pdf,dataBoth,frBoth,
                        categoryTitle,energyStr,lumiDict[energyStr],
                        caption2=r"DY+t\bar{t} MC",
                        showStackDatasets=[dataTT,dataDY],
                        showStackDatasetTitles=[legendEntries['ttbar'],legendEntries["DYJetsToLL"]],
                        showStackDatasetColors=[colors['ttbar'],colors["DYJetsToLL"]],
                        errorsToDraw=errorsToDraw,
                        canvas=canvas
                        )
  rmpBoth.draw("mcFit_"+category+"_Both")
  
  #####################################################################
  # Vary ttbar up and down
  
  # ttbar varied up
  
  factorTT = 5.
  print "vary ttbar up by: ",factorTT
  weightTT = lumiDict[energyStr]*1000.*xsec["ttbar_"+energyStr]/nEventsMap["ttbar_"+energyStr]*factorTT
  
  weightVar.setVal(weightTT)
  dataTTUp = root.RooDataSet("dataTTUp","TT Data Scaled Up",root.RooArgSet(dimuonMass,weightVar),root.RooFit.Import(dataNoWeightTT),root.RooFit.WeightVar(weightVar))
  weightVar.setVal(weightDY)
  dataBothTTUp = root.RooDataSet("dataBothTTUp","DY+TTUp Data",dataDY,root.RooArgSet(dimuonMass,weightVar),"","weight")
  weightVar.setVal(weightTT)
  dataBothTTUp.append(dataTTUp)
  
  dataTTUp.Print()
  dataBothTTUp.Print()
  
  frBothTTUp = pdf.fitTo(dataBothTTUp,root.RooFit.Save(),PRINTLEVEL,root.RooFit.SumW2Error(useWeight2Error))
#  frBothTTUp.Print()
  rmpBothTTUp = RooModelPlotter(dimuonMass,pdf,dataBothTTUp,frBothTTUp,
                        categoryTitle,energyStr,lumiDict[energyStr],
                        caption2=r"DY+t\bar{t} MC, t\bar{t} Scaled Up",
                        showStackDatasets=[dataTTUp,dataDY],
                        showStackDatasetTitles=[legendEntries['ttbar']+" #times {0}".format(factorTT),legendEntries["DYJetsToLL"]],
                        showStackDatasetColors=[colors['ttbar'],colors["DYJetsToLL"]],
                        errorsToDraw=errorsToDraw,
                        canvas=canvas
                        )
  rmpBothTTUp.draw("mcFit_"+category+"_BothTTUp")
  
  
  # ttbar varied Down
  
  factorTT = 0.2
  print "vary ttbar down by: ",factorTT
  weightTT = lumiDict[energyStr]*1000.*xsec["ttbar_"+energyStr]/nEventsMap["ttbar_"+energyStr]*factorTT
  
  weightVar.setVal(weightTT)
  dataTTDown = root.RooDataSet("dataTTDown","TT Data Scaled Down",root.RooArgSet(dimuonMass,weightVar),root.RooFit.Import(dataNoWeightTT),root.RooFit.WeightVar(weightVar))
  weightVar.setVal(weightDY)
  dataBothTTDown = root.RooDataSet("dataBothTTDown","DY+TTDown Data",dataDY,root.RooArgSet(dimuonMass,weightVar),"","weight")
  weightVar.setVal(weightTT)
  dataBothTTDown.append(dataTTDown)
  
  dataTTDown.Print()
  dataBothTTDown.Print()
  
  frBothTTDown = pdf.fitTo(dataBothTTDown,root.RooFit.Save(),PRINTLEVEL,root.RooFit.SumW2Error(useWeight2Error))
#  frBothTTDown.Print()
  rmpBothTTDown = RooModelPlotter(dimuonMass,pdf,dataBothTTDown,frBothTTDown,
                        categoryTitle,energyStr,lumiDict[energyStr],
                        caption2=r"DY+t\bar{t} MC, t\bar{t} Scaled Down",
                        showStackDatasets=[dataTTDown,dataDY],
                        showStackDatasetTitles=[legendEntries['ttbar']+" #times {0}".format(factorTT),legendEntries["DYJetsToLL"]],
                        showStackDatasetColors=[colors['ttbar'],colors["DYJetsToLL"]],
                        errorsToDraw=errorsToDraw,
                        canvas=canvas
                        )
  rmpBothTTDown.draw("mcFit_"+category+"_BothTTDown")
