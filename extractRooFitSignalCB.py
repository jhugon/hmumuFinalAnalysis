#! /usr/bin/env python

from ROOT import gSystem

import sys
import os
import re
import math
from ROOT import *
gSystem.Load('libRooFit')
import ROOT as root

#from helpers import *

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT


# ------------------------------------------------------------------------------

def Usage():
   print 'Wrong syntax: \n'
   print '   ./extractRooFitSignalDG.py [process] [mass] [benergy (7TeV,8TeV)]\n'
   print 'first is included, last IS NOT included'
   sys.exit()
    
# ------------------------------------------------------------------------------

# ====================   
# set useful variables
if (len(sys.argv) < 4):
   Usage()

# default values
process = sys.argv[1]
mass    = sys.argv[2]
benergy = sys.argv[3]

#inputfile  = "/afs/cern.ch/work/d/digiovan/HToMM/officialMC/submit%s/root/" % benergy
inputfile  = "/afs/cern.ch/work/d/digiovan/HToMM/baselinePP/submit%s/root/" % benergy
inputfile += "%sHmumu%s_%s.root" % (process,mass,benergy)

if (process == 'GluGlu'):
   process = 'gg'
if (process == 'VBF'):
   process = 'vbf'
   
outputfile = 'fitresults/fitresults_%sHmumu%s_%s' % (process,mass,benergy)
savefile   = 'fitresults/png/fitresults_%sHmumu%s_%s' % (process,mass,benergy)

print 'PROCESS: %s' % process
print 'MASS: %s' % mass
print 'BENERGY: %s' % benergy
print 'inputfile: %s' % inputfile
print 'outputfile: %s' % outputfile
print 'savefile: %s' % savefile
print '----------------------------------------------------------------------\n'
# ====================


# define the dimuon mass variable
#minMass = 110.
#maxMass = 135.
minMass = float(mass.replace("p","."))-10.0
maxMass = float(mass.replace("p","."))+5.0
mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)

#####################################################################
# define the Crystal Ball
mean   = root.RooRealVar("Mean","Mean",float(mass),float(mass)-5,float(mass)+5)

width  = root.RooRealVar("Width","Width",1.5,0.5,5.0)

Alpha = root.RooRealVar("Alpha","Alpha",1.8,-10.0,10.0)
n     = root.RooRealVar("n","n",1.0,-50.0,50.0)

pdfMmumuCB = root.RooCBShape ("pdf","pdf",mMuMu,mean,width,Alpha,n)

#####################################################################


def doFit(f,baseName,canvas,outputfile,saveName):

  mDiMu = f.Get( baseName + "/mDiMu")

  #mDiMu.Rebin(2)
  
  mMuMuRooDataHist = root.RooDataHist("template","template",root.RooArgList(mMuMu),mDiMu)

  pdfMmumu = pdfMmumuCB
    
  pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.SumW2Error(False),PRINTLEVEL)
  
  canvas.cd()
  #canvas.SetLogy()
  
  plotMmumu = mMuMu.frame(minMass,maxMass)
  
  mMuMuRooDataHist.plotOn(plotMmumu)
  
  pdfMmumu .plotOn(plotMmumu,root.RooFit.LineColor(root.kAzure+2))

  plotMmumu.SetTitle(baseName)
  if ('hists' in baseName):
     plotMmumu.SetTitle(baseName[5:])
     
  plotMmumu.Draw()

##  chi2ondf =  plotMmumu.chiSquare()
##  tlatex = root.TLatex()
##  tlatex.SetNDC()
##  tlatex.SetTextFont(root.gStyle.GetLabelFont())
##  tlatex.SetTextSize(0.05)
##  tlatex.SetTextAlign(22)
##  tlatex.DrawLatex(0.75,0.85,"#chi^{{2}}/NDF = {0:.2f}".format(chi2ondf))

  canvas.Update()
  if saveName != "":
     canvas.SaveAs(saveName)
 
  outfile = open(outputfile+'_CB','w')
          
  print ''
  print ''
  print 'Crystal Ball'
  for i in [mean,width,Alpha,n]:
     print("{0}: {1:.10g}".format(i.GetName(),i.getVal()))
      
  outfile.write('#process mass mean err width err alpha err n err\n')
  outfile.write('%s %s %s %s %s %s %s %s %s %s\n'
                % (process,mass,
                     mean.getVal(),   mean.getError(), 
                    width.getVal(),  width.getError(),
                    Alpha.getVal(),  Alpha.getError(),
                        n.getVal(),      n.getError()
                   )
                )

  
#####################################################################

f = root.TFile(inputfile)

# baseline
#baseNamesGG = ["IncPreselPtG10BB",
#               "IncPreselPtG10BO",
#               "IncPreselPtG10BE",
#               "IncPreselPtG10OO",
#               "IncPreselPtG10OE",
#               "IncPreselPtG10EE"
#               ]
#
#baseNamesVBF = ["VBFBDTCut",
#                "histsVBFDeJJG3p5MJJG550pTmissL100",
#                "histsVBFDeJJG3p4MJJG500pTmissL25"
#                ]
            
# baseline++

baseNamesGG = ["Jets01PassPtG10BB",
               "Jets01PassPtG10BO",
               "Jets01PassPtG10BE",
               "Jets01PassPtG10OO",
               "Jets01PassPtG10OE",
               "Jets01PassPtG10EE",

               "Jets01FailPtG10BB",
               "Jets01FailPtG10BO",
               "Jets01FailPtG10BE",
               "Jets01FailPtG10OO",
               "Jets01FailPtG10OE",
               "Jets01FailPtG10EE",               
               ]

baseNamesVBF = ["Jet2CutsVBFPass",  
                "Jet2CutsGFPass",   
                "Jet2CutsFailVBFGF"
                ]

baseNames = None
if ("GluGluHmumu" in f.GetName()):
  baseNames = baseNamesGG
if ("VBFHmumu" in f.GetName()):
  baseNames = baseNamesVBF

if (baseNames == None):
  print "ERROR: The input file is not ggH or vbfH"
  sys.exit()

  
canvases = []

iStep = 0
for baseName in baseNames:
  canvas = root.TCanvas("canvas_%s" % baseName,"",iStep,iStep,800,600)
  canvases.append( canvas )
  iStep += 100
  
for id in range(0,len(baseNames)):
   if ('hists' in baseNames[id]):
      doFit(f,baseNames[id],canvases[id],outputfile+'_'+baseNames[id][5:],savefile+'_'+baseNames[id][5:]+'_CB.png')
   else:
      doFit(f,baseNames[id],canvases[id],outputfile+'_'+baseNames[id],savefile+'_'+baseNames[id]+'_CB.png')
