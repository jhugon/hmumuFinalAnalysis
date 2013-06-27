#! /usr/bin/env python

from ROOT import gSystem

import sys
import os
import re
import math
import time
from ROOT import *
gSystem.Load('libRooFit')
import ROOT as root

from helpers import *

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT


# ------------------------------------------------------------------------------

def Usage():
   print 'Wrong syntax: \n'
   print '   ./muonMomScaleSyst.py [process] [benergy (7TeV,8TeV)]\n'
   print 'first is included, last IS NOT included'
   sys.exit()
    
# ------------------------------------------------------------------------------

# ====================   
# set useful variables
if (len(sys.argv) < 3):
   Usage()

# default values
process = sys.argv[1]
benergy = sys.argv[2]

print ' '
print 'PROCESS: %s' % process
print 'BENERGY: %s' % benergy
print '----------------------------------------------------------------------\n'


# ====================
def doSystVariation(tree,baseName,canvas,process,mass,benergy):

   #minMass =  70.
   #maxMass = 170.
   
   minMass = float(mass.replace("p","."))-15.0
   maxMass = float(mass.replace("p","."))+10.0

   # 
   category  = baseName[0]
   cutString = baseName[1]
   fullCutString = treeCut(category,cutString,eventWeights=False,muonRequirements=True)
   
   #print "analyzing category \"%s\"" % category
   #print "applying cut \"%s\"" % fullCutString

   dimuonMass          = root.RooRealVar("dimuonMass",         "M(#mu#mu) [GeV/c^{2}]",minMass,maxMass)
   dimuonMassResSFUp   = root.RooRealVar("dimuonMassResSFUp",  "M(#mu#mu) [GeV/c^{2}]",minMass,maxMass)
   dimuonMassResSFDown = root.RooRealVar("dimuonMassResSFDown","M(#mu#mu) [GeV/c^{2}]",minMass,maxMass)

   # D A T A S E T S
   mMuMuRooDataSets = []

   # default
   rooDataSet = root.RooDataSet("roodataset_%s" % category, "",
                                tree, root.RooArgSet(dimuonMass))

   mMuMuRooDataSets += [[rooDataSet, dimuonMass, root.kBlack, "Baseline"]]

   # muon momentum scale resolution up
   rooDataSetResSFUp = root.RooDataSet("roodatasetresup_%s" % category, "",
                                       tree, root.RooArgSet(dimuonMassResSFUp))

   mMuMuRooDataSets += [[rooDataSetResSFUp, dimuonMassResSFUp, root.kRed, "ResUp"]]


   # muon momentum scale resolution down
   rooDataSetResSFDown = root.RooDataSet("roodatasetresdown_%s" % category, "",
                                         tree, root.RooArgSet(dimuonMassResSFDown))

   mMuMuRooDataSets += [[rooDataSetResSFDown, dimuonMassResSFDown, root.kAzure+2, "ResDown"]]



   # F I T   and  D R A W
   canvas.Clear()
   canvas.cd()

   pars = []
   # F I T
   for dataset in mMuMuRooDataSets:

      print dataset
      mMuMu     = dataset[1]
      plotMmumu = mMuMu.frame(minMass,maxMass)
      plotMmumu.SetTitle( ("%s, " + category + " m=%s GeV/c^{2} at %s") % (process,mass,benergy) )

      # plot the dataset
      dataset[0].plotOn( plotMmumu,
                         root.RooFit.Binning(50),
                         root.RooFit.MarkerColor(dataset[2]),
                         root.RooFit.LineColor(dataset[2])
                       )
      
      
      # F I T   F U N C T I O N
      # define the Double Gaussian

      meanG1 = root.RooRealVar("MeanG1","MeanG1", float(mass), float(mass)-10, float(mass)+3 )
      meanG2 = root.RooRealVar("MeanG2","MeanG2", float(mass), float(mass)-3,  float(mass)+3 )
      
      widthG1 = root.RooRealVar("WidthG1","WidthG1", 7.049779267, 0.0,50.0)
      widthG2 = root.RooRealVar("WidthG2","WidthG2", 1.830513636, 0.0, 4.0)
      
      mixGG = root.RooRealVar("mixGG","mixGG", 0.1140210709, 0.0,1.0)
      
      gaus1 = root.RooGaussian("gaus1","gaus1",mMuMu,meanG1,widthG1)
      gaus2 = root.RooGaussian("gaus2","gaus2",mMuMu,meanG2,widthG2)
      
      pdfMmumuGG = root.RooAddPdf("pdfMmumuGG","pdfMmumuGG",gaus1,gaus2,mixGG)


      # fit the dataset
      pdfMmumuGG.fitTo(dataset[0],root.RooFit.SumW2Error(False),PRINTLEVEL)
      pdfMmumuGG.plotOn(plotMmumu,root.RooFit.LineColor(dataset[2]))


      mean_narrow      = meanG2.getVal()
      err_mean_narrow  = meanG2.getError()
      width_narrow     = widthG2.getVal()
      err_width_narrow = widthG2.getError()

      mean_wide        = meanG1.getVal()
      err_mean_wide    = meanG1.getError()
      width_wide       = widthG1.getVal()
      err_width_wide   = widthG1.getError()

      mixing           = mixGG.getVal()
      err_mixing       = mixGG.getError()

      if (widthG1.getVal()<widthG2.getVal()):

         mean_narrow      = meanG1.getVal()
         err_mean_narrow  = meanG1.getError()
         width_narrow     = widthG1.getVal()
         err_width_narrow = widthG1.getError()
         
         mean_wide        = meanG2.getVal()
         err_mean_wide    = meanG2.getError()
         width_wide       = widthG2.getVal()
         err_width_wide   = widthG2.getError()
         
         mixing           = (1-mixGG.getVal())
         err_mixing       = mixGG.getError()

      print "AAAAHOOOOOO: ", [mean_narrow,mean_wide,width_narrow,width_wide,mixing]
      pars += [ [mean_narrow,mean_wide,width_narrow,width_wide,mixing] ]

      plotMmumu.Draw()
      canvas.Update()
      canvas.SaveAs("systematics/muonMomScaleSyst/" + process + "_" + category + "_" + mass + "_" + dataset[3] + ".png")

   # return relative yield changes
   percUp   = abs( (pars[0][2] - pars[1][2]) / pars[0][2] )
   percDown = abs( (pars[0][2] - pars[2][2]) / pars[0][2] )

   print float( "%.6f" % max(percUp, percDown) )
   return float( "%.6f" % max(percUp, percDown) )



#####################################################################


# baseline++
jet2PtCuts = " jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40."
jet01PtCuts = " !(jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40.)"

baseNames = []

baseNames += [["Jets01PassPtG10",    "dimuonPt>10. && " + jet01PtCuts]]
baseNames += [["Jets01FailPtG10","  !(dimuonPt>10.) && "+ jet01PtCuts]]
baseNames += [["Jet2CutsVBFPass",    "deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts]]
baseNames += [["Jet2CutsGFPass",   "!(deltaEtaJets>3.5 && dijetMass>650.) &&  (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
baseNames += [["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]


masses = ['115', '125', '135', '150']
#masses = ['125']

hSyst = root.TH2F("hSyst", "%s at %s" % (process, benergy),
                  len(masses),0,len(masses),
                  len(baseNames),0,len(baseNames) )

for id in range(0,len(masses)):
   hSyst.GetXaxis().SetBinLabel(id+1,masses[id])

for id in range(0,len(baseNames)):
   hSyst.GetYaxis().SetBinLabel(id+1,baseNames[id][0])


folder  = "/afs/cern.ch/work/d/digiovan/HToMM/baselinePP/submit%s/root/" % benergy
#folder  = "/afs/cern.ch/work/d/digiovan/code/CMSSW_6_1_1/results/scripts/htomm/justin/development/hmumuAnalysisV00-01-10/ptM15GeV/"
for massId in range(0,len(masses)):

   mass = masses[massId]
   inputfile = folder + "%sHmumu%s_%s.root" % (process,mass,benergy)

   print 'Opening file %s' % inputfile
   file = root.TFile(inputfile)

   canvases = []

   iStep = 0
   for baseName in baseNames:
      canvas = root.TCanvas("canvas_%s_m%s" % (baseName[0],mass),"",iStep,iStep,800,600)
      canvases.append( canvas )
      iStep += 100

   for bId in range(0,len(baseNames)):

      treename = "outtree"+str(baseNames[bId][0])
      print "Getting TTree %s" % treename
      tree = file.Get( treename )
      tree.SetCacheSize(10000000);
      tree.AddBranchToCache("*");

      #doSystVariation(tree,baseNames[bId],canvases[bId])
      syst = doSystVariation(tree,baseNames[bId],canvases[bId],process,mass,benergy)
      hSyst.SetBinContent(massId+1,bId+1,syst)



canvas = root.TCanvas("canvas","",300,300,800,600)
canvas.SetMargin(0.19,0.05,0.13,0.08)
canvas.cd()

hSyst.Draw("COLZ TEXT")

canvas.Update()
canvas.SaveAs("systematics/muonMomScaleSyst/syst_muonMomRes_" + process + "_" + benergy +".png")
canvas.SaveAs("systematics/muonMomScaleSyst/syst_muonMomRes_" + process + "_" + benergy +".pdf")


