#! /usr/bin/env python

from ROOT import gSystem

import sys
import os
import re
import math
import commands
from ROOT import *
gSystem.Load('libRooFit')
import ROOT as root

from helpers import *

# ------------------------------------------------------------------------------

def Usage():
   print 'Wrong syntax: \n'
   print '   ./plotRooFitSignal.py [process] [benergy (7TeV,8TeV)] [fit (DG,CBG)] [cat (BB,BO,...,EE,VBF)]\n'
   print 'first is included, last IS NOT included'
   sys.exit()
    
# ------------------------------------------------------------------------------

# ====================   
# set useful variables
if (len(sys.argv) < 5):
   Usage()

# default values
process = sys.argv[1]
benergy = sys.argv[2]
fit     = sys.argv[3]
cat     = sys.argv[4]

folder  = "fitresults/"

print 'PROCESS: %s' % process
print 'BENERGY: %s' % benergy
print '----------------------------------------------------------------------\n'
# ====================


#####################################################################


baseName = None
title = None

if (process == 'VBF' and 'VBF' not in cat):
   print "Mismatch -> exiting...\n"
   sys.exit()
   
if (process == 'GluGlu' and 'VBF' in cat):
   print "Mismatch -> exiting...\n"
   sys.exit()


if (cat == "BB"):
   baseName = "IncPreselPtG10BB"
   title = "gg, Barrel - Barrel (%s)" % benergy
if (cat == "BO"):
   baseName = "IncPreselPtG10BO"
   title = "gg, Barrel - Overlap (%s)" % benergy
if (cat == "BE"):
   baseName = "IncPreselPtG10BE"
   title = "gg, Barrel - Endcap (%s)" % benergy
if (cat == "OO"):
   baseName = "IncPreselPtG10OO"
   title = "gg, Overlap - Overlap (%s)" % benergy
if (cat == "OE"):
   baseName = "IncPreselPtG10OE"
   title = "gg, Overlap - Endcap (%s)" % benergy
if (cat == "EE"):
   baseName = "IncPreselPtG10EE"
   title = "gg, Endcap - Endcap (%s)" % benergy
      
if (cat == "VBF"):
   baseName = "VBFBDTCut"
   title = "VBF (%s)" % benergy
if (cat == "VBFCiCLoose"):
   baseName = "histsVBFDeJJG3p5MJJG550pTmissL100"
   title = "VBF CiC Loose (%s)" % benergy
if (cat == "VBFCiCTight"):
   baseName = "histsVBFDeJJG3p4MJJG500pTmissL25"
   title = "VBF CiC Tight (%s)" % benergy

if (baseName == None):
   print "No category defined -> exiting...\n"
   sys.exit()
   

varNames = ['meanG1','widthG1','meanG2','widthG2','mixGG']

listOfFiles = commands.getoutput('ls %s | grep %s | grep %s | grep %s | grep fit' % (folder,baseName,benergy,fit) ).split('\n')
print listOfFiles

  
canvases = []
graphs = []

iStep = 0
for varName in varNames:
  canvas = root.TCanvas("canvas_%s" % varName,"",iStep,iStep,800,600)
  canvases.append( canvas )
  iStep += 100

  graph = root.TGraphErrors(len(listOfFiles))
  graphs.append ( graph )


index = 0
for file in listOfFiles:
   input = open(folder+"/"+file)

   for line in input:
      if ('#' in line):
         continue
      else:
         line = line.split('\n')[0]
         line = line.split(' ')
         #print line

         process = line[0]
         mass    = line[1]

         for id in range(0,len(varNames)):
            parValue = line[2+(2*id)]
            parError = line[2+(2*id)+1]
            
            #print index, mass, 
            graphs[id].SetPoint(int(index),float(mass),float(parValue))
            graphs[id].SetPointError(int(index),float(0),float(parError))
            
         index+=1


pathToSave = 'fitresults/'

for id in range(0,len(varNames)):
   canvases[id].cd()
   graphs[id].GetXaxis().SetTitle("Mass(#mu#mu)")
   graphs[id].GetYaxis().SetTitle(varNames[id])
   graphs[id].SetTitle(title)
   graphs[id].Draw("AP")

   fitFunc = root.TF1("fitFunc","pol1",120.,140.)
   graphs[id].Fit("fitFunc")

   legend = root.TLegend(0.56,0.76,0.66,0.87)
   legend.SetFillColor(0)
   legend.SetBorderSize(0)
   legend.SetTextSize(0.042)
   legtitle = "pol1, #chi^{{2}}/ndf= {0:.2f}/{1}".format(fitFunc.GetChisquare(),
                                                         fitFunc.GetNDF())
   legend.AddEntry(fitFunc,
                   legtitle,
                   "l")
   
   legend.Draw("same")
   
   canvases[id].Update()

   saveName = pathToSave + 'png/' + varNames[id] + "_" + process + '_' + benergy + '_' +  cat + '.png'
   canvases[id].SaveAs(saveName)

   # save the output of the fit
   outputfile = pathToSave + 'extrapolation_' + varNames[id] + "_" + process + '_' + benergy + '_' +  cat 
   outfile = open(outputfile+'_DG.txt','w')

   outfile.write('#par0 err_par0 par1 err_par1 chisquare ndf\n')  
   outfile.write('{0:.10g} {1:.10g} {2:.10g} {3:.10g} {4:.4g} {5}\n'.format(fitFunc.GetParameter(0),
                                                                            fitFunc.GetParError(0),
                                                                            fitFunc.GetParameter(1),
                                                                            fitFunc.GetParError(1),
                                                                            fitFunc.GetChisquare(),
                                                                            fitFunc.GetNDF()
                                                                            )
                 ) 


