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

inputfile  = "/afs/cern.ch/work/d/digiovan/HToMM/officialMC/submit%s/root/" % benergy
inputfile += "%sHmumu%s_%s.root" % (process,mass,benergy)

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
minMass = float(mass)-10.0
maxMass = float(mass)+10.0
mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)

#####################################################################
# define the Double Gaussian
meanG1 = root.RooRealVar("MeanG1","MeanG1", float(mass)-4.0, 100.,170.)
meanG2 = root.RooRealVar("MeanG2","MeanG2", float(mass), 100.,170.)
 
widthG1 = root.RooRealVar("WidthG1","WidthG1", 5.049779267, 0.1,20.0)
widthG2 = root.RooRealVar("WidthG2","WidthG2", 1.830513636, 0.1,5.0)
 
mixGG = root.RooRealVar("mixGG","mixGG", 0.1140210709, 0.0,1.0)

gaus1 = root.RooGaussian("gaus1","gaus1",mMuMu,meanG1,widthG1)
gaus2 = root.RooGaussian("gaus2","gaus2",mMuMu,meanG2,widthG2)
 
pdfMmumuGG = root.RooAddPdf("pdfMmumuGG","pdfMmumuGG",gaus1,gaus2,mixGG)


# define the Single Gaussian for the EE category
meanSG  = root.RooRealVar("MeanSG", "MeanSG", float(mass),110.,140.)
widthSG = root.RooRealVar("WidthSG","WidthSG", 1.839226054,0.1,20.0)
 
pdfMmumuSG = root.RooGaussian("sgaus","sgaus",mMuMu,meanSG,widthSG)


#####################################################################


def doFit(f,baseName,canvas,outputfile,saveName):

  mDiMu = f.Get( baseName + "/mDiMu")

  #mDiMu.Rebin(2)
  
  mMuMuRooDataHist = root.RooDataHist("template","template",root.RooArgList(mMuMu),mDiMu)

  pdfMmumu = None
#  if (baseName == "IncPreselPtG10EE"):
#    pdfMmumu = pdfMmumuSG
#  else:
  pdfMmumu = pdfMmumuGG
    
    
  pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.SumW2Error(False),PRINTLEVEL)
  
  canvas.cd()
  #canvas.SetLogy()
  
  plotMmumu = mMuMu.frame(minMass,maxMass)
  
  mMuMuRooDataHist.plotOn(plotMmumu)
  
  pdfMmumu .plotOn(plotMmumu,root.RooFit.LineColor(root.kAzure+2))
  
  plotMmumu.SetTitle(baseName)
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
 
  outfile = open(outputfile+'_DG','w')
          
  print ''
  print ''
#  if (baseName == "IncPreselPtG10EE"):
#    print 'Single Gaussian'
#    for i in [meanSG,widthSG]:
#      print("{0}: {1:.10g}".format(i.GetName(),i.getVal()))
#      
#    outfile.write('#process mass meanSG err widthSG err\n')
#    outfile.write('%s %s %s %s %s %s\n' % (process,mass,
#                                           meanSG.getVal() , meanSG.getError(), 
#                                           widthSG.getVal(), widthSG.getError()
#                                           )
#                  )
#  else:
  print 'Double Gaussian'
  for i in [meanG1,meanG2,widthG1,widthG2,mixGG]:
     print("{0}: {1:.10g}".format(i.GetName(),i.getVal()))

  outfile.write('#process mass meanG1 err widthG1 err meanG2 err widthG2 err mixGG err\n')  
  outfile.write('%s %s %s %s %s %s %s %s %s %s %s %s\n'
                % (process,mass,
                   meanG1.getVal(),  meanG1.getError(), 
                   widthG1.getVal(), widthG1.getError(),
                   meanG2.getVal(),  meanG2.getError(),
                   widthG2.getVal(), widthG2.getError(),
                   mixGG.getVal(),   mixGG.getError()
                   )
                ) 

  
#####################################################################

f = root.TFile(inputfile)

baseNamesGG = ["IncPreselPtG10BB",
               "IncPreselPtG10BO",
               "IncPreselPtG10BE",
               "IncPreselPtG10OO",
               "IncPreselPtG10OE",
               "IncPreselPtG10EE"
               ]

baseNamesVBF = ["VBFBDTCut",
                "histsVBFDeJJG3p5MJJG550pTmissL100",
                "histsVBFDeJJG3p4MJJG500pTmissL25"
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
  doFit(f,baseNames[id],canvases[id],outputfile+'_'+baseNames[id],savefile+'_'+baseNames[id]+'.png')
