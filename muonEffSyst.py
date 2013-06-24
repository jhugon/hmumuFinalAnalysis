#! /usr/bin/env python

from ROOT import gSystem

import sys
import os
import re
import math
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
   print '   ./muonEffSyst.py [process (GluGlu,VBF)] [benergy (7TeV,8TeV)]\n'
   sys.exit()
    
# ------------------------------------------------------------------------------

# ====================   
# set useful variables
if (len(sys.argv) < 3):
   Usage()

# default values
process = sys.argv[1]
benergy = sys.argv[2]

print 'PROCESS: %s' % process
print 'BENERGY: %s' % benergy
print '----------------------------------------------------------------------\n'
# ====================

def doSystVariation(tree,baseName,histos,canvas):

   category  = baseName[0]
   cutString = baseName[1]

   fullCutString = treeCut(category,cutString,eventWeights=False,muonRequirements=True)

   #print "analyzing category \"%s\"" % category
   #print "applying cut \"%s\"" % fullCutString

   hBase = histos[0]
   hUp   = histos[1]
   hDown = histos[2]
   
   drawString = "dimuonMass >> {0}".format(hBase.GetName())
   tree.Draw(drawString, "("+fullCutString+")*puWeight" )

   drawStringUp = "dimuonMass >> {0}".format(hUp.GetName())
   tree.Draw(drawStringUp, "("+fullCutString+")*puWeightMuonEffUp" )

   drawStringDown = "dimuonMass >> {0}".format(hDown.GetName())
   tree.Draw(drawStringDown, "("+fullCutString+")*puWeightMuonEffDown" )

   # Draw
   canvas.cd()

   hBase.SetTitle(category)
   hBase.GetXaxis().SetTitle("M(#mu#mu) [GeV/c^{2}]")
   hBase.Draw()
   hUp  .Draw("same")
   hDown.Draw("same")

   canvas.Update()

   # return relative yield changes
   percUp   = abs(   (hUp.Integral() - hBase.Integral()) / hBase.Integral())
   percDown = abs( (hDown.Integral() - hBase.Integral()) / hBase.Integral())

   return float( "%.6f" % max(percUp, percDown) )

#####################################################################

           
# baseline++

jet2PtCuts = " && jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40."
jet01PtCuts = " && !(jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40.)"

baseNames = []

baseNames += [["Jets01PassPtG10",    "dimuonPt>10." + jet01PtCuts]]
baseNames += [["Jets01FailPtG10","  !(dimuonPt>10.)"+ jet01PtCuts]]
baseNames += [["Jet2CutsVBFPass",    "deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts]]
baseNames += [["Jet2CutsGFPass",   "!(deltaEtaJets>3.5 && dijetMass>650.) &&  (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
baseNames += [["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]


masses = ['115', '125', '135', '150']
hSyst = root.TH2F("hSyst", "%s at %s" % (process, benergy),
                  len(masses),0,len(masses),
                  len(baseNames),0,len(baseNames) )

for id in range(0,len(masses)):
   hSyst.GetXaxis().SetBinLabel(id+1,masses[id])

for id in range(0,len(baseNames)):
   hSyst.GetYaxis().SetBinLabel(id+1,baseNames[id][0])





folder  = "/afs/cern.ch/work/d/digiovan/HToMM/baselinePP/submit%s/root/" % benergy
for massId in range(0,len(masses)):

   mass = masses[massId]
   inputfile = folder + "%sHmumu%s_%s.root" % (process,mass,benergy)

   file = root.TFile(inputfile)

   tree = file.Get("outtree")
   tree.SetCacheSize(10000000);
   tree.AddBranchToCache("*");


   canvases = []
   histos   = []

   iStep = 0
   for baseName in baseNames:
      canvas = root.TCanvas("canvas_%s_m%s" % (baseName[0],mass),"",iStep,iStep,800,600)
      canvases.append( canvas )
      iStep += 100

      mDiMu     = root.TH1F("mDiMu_%s_m%s"     % (baseName[0],mass), "",60,110,170)
      mDiMuUp   = root.TH1F("mDiMuUp_%s_m%s"   % (baseName[0],mass), "",60,110,170)
      mDiMuDown = root.TH1F("mDiMuDown_%s_m%s" % (baseName[0],mass), "",60,110,170)

      mDiMu     . SetLineColor(root.kBlack)
      mDiMuUp   . SetLineColor(root.kRed)
      mDiMuDown . SetLineColor(root.kBlue)

      histos += [ [mDiMu,mDiMuUp,mDiMuDown] ]


   for bId in range(0,len(baseNames)):

      syst = doSystVariation(tree,baseNames[bId],histos[bId],canvases[bId])
      hSyst.SetBinContent(massId+1,bId+1,syst)



canvas = root.TCanvas("canvas","",300,300,800,600)
canvas.SetMargin(0.19,0.05,0.13,0.08)
canvas.cd()
hSyst.Draw("COLZ TEXT")
canvas.Update()
