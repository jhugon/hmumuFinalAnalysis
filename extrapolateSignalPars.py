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


if (cat == "Jets01PassPtG10BB"):
   baseName = "Jets01PassPtG10BB"
   title = "Jet 0+1 cat p_{T}(#mu#mu)>10 GeV/c, Barrel - Barrel (%s)" % benergy
if (cat == "Jets01PassPtG10BO"):
   baseName = "Jets01PassPtG10BO"
   title = "Jet 0+1 cat p_{T}(#mu#mu)>10 GeV/c, Barrel - Overlap (%s)" % benergy
if (cat == "Jets01PassPtG10BE"):
   baseName = "Jets01PassPtG10BE"
   title = "Jet 0+1 cat p_{T}(#mu#mu)>10 GeV/c, Barrel - Endcap (%s)" % benergy
if (cat == "Jets01PassPtG10OO"):
   baseName = "Jets01PassPtG10OO"
   title = "Jet 0+1 cat p_{T}(#mu#mu)>10 GeV/c, Overlap - Overlap (%s)" % benergy
if (cat == "Jets01PassPtG10OE"):
   baseName = "Jets01PassPtG10OE"
   title = "Jet 0+1 cat p_{T}(#mu#mu)>10 GeV/c, Overlap - Endcap (%s)" % benergy
if (cat == "Jets01PassPtG10EE"):
   baseName = "Jets01PassPtG10EE"
   title = "Jet 0+1 cat p_{T}(#mu#mu)>10 GeV/c, Endcap - Endcap (%s)" % benergy
if (cat == "Jets01PassPtG10CC"):
   baseName = "Jets01PassPtG10CC"
   title = "Jet 0+1 cat p_{T}(#mu#mu)>10 GeV/c, BE+OO (%s)" % benergy
if (cat == "Jets01PassPtG10FF"):
   baseName = "Jets01PassPtG10FF"
   title = "Jet 0+1 cat p_{T}(#mu#mu)>10 GeV/c, OE+EE (%s)" % benergy
      
if (cat == "Jets01FailPtG10BB"):
   baseName = "Jets01FailPtG10BB"
   title = "Jet 0+1 cat p_{T}(#mu#mu)<10 GeV/c, Barrel - Barrel (%s)" % benergy
if (cat == "Jets01FailPtG10BO"):
   baseName = "Jets01FailPtG10BO"
   title = "Jet 0+1 cat p_{T}(#mu#mu)<10 GeV/c, Barrel - Overlap (%s)" % benergy
if (cat == "Jets01FailPtG10BE"):
   baseName = "Jets01FailPtG10BE"
   title = "Jet 0+1 cat p_{T}(#mu#mu)<10 GeV/c, Barrel - Endcap (%s)" % benergy
if (cat == "Jets01FailPtG10OO"):
   baseName = "Jets01FailPtG10OO"
   title = "Jet 0+1 cat p_{T}(#mu#mu)<10 GeV/c, Overlap - Overlap (%s)" % benergy
if (cat == "Jets01FailPtG10OE"):
   baseName = "Jets01FailPtG10OE"
   title = "Jet 0+1 cat p_{T}(#mu#mu)<10 GeV/c, Overlap - Endcap (%s)" % benergy
if (cat == "Jets01FailPtG10EE"):
   baseName = "Jets01FailPtG10EE"
   title = "Jet 0+1 cat p_{T}(#mu#mu)<10 GeV/c, Endcap - Endcap (%s)" % benergy
if (cat == "Jets01FailPtG10CC"):
   baseName = "Jets01FailPtG10CC"
   title = "Jet 0+1 cat p_{T}(#mu#mu)<10 GeV/c, BE+OO (%s)" % benergy
if (cat == "Jets01FailPtG10FF"):
   baseName = "Jets01FailPtG10FF"
   title = "Jet 0+1 cat p_{T}(#mu#mu)<10 GeV/c, OE+EE (%s)" % benergy


if (cat == "Jet2CutsVBFPass"):
   baseName = "Jet2CutsVBFPass"
   title = ">=2 Jets, VBF Pass (%s)" % benergy
if (cat == "Jet2CutsGFPass"):
   baseName = "Jet2CutsGFPass"
   title = ">=2 Jets, VBF Fail, GG+VH Pass (%s)" % benergy
if (cat == "Jet2CutsFailVBFGF"):
   baseName = "Jet2CutsFailVBFGF"
   title = ">=2 Jets, VBF Fail, GG+VH Fail (%s)" % benergy

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

         #print index
         print line[2+(2*(len(varNames)-1))]
         if ( float(line[2+(2*(len(varNames)-1))]) > 0.9 or float(line[2+(2*(len(varNames)-1))]) < 0.05):
            continue
         
         for id in range(0,len(varNames)):
            parValue = line[2+(2*id)]
            parError = line[2+(2*id)+1]
            
            #print  parValue,parError
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

   fitFunc = root.TF1("fitFunc","pol1",110.,160.)
   #fitFunc = root.TF1("fitFunc","[0]*x",110.,160.)
   #fitFunc.FixParameter(0,0)
   #fitFunc.SetParLimits(0,1,0)
   #fitFunc.SetParameter(1,0.1)

   # The "Fail" categories unfortunately need some hacking...
   #if ("Fail" in cat and id==len(varNames)-1):
   #
   #   fitFunc.SetParameter(0,0.1)
   #   if ("BE" in cat or "OO" in cat or "OE" in cat or "FF" in cat):
   #      fitFunc.SetParLimits(0,3,2)
   #
   #   if ("EE" in cat):
   #      fitFunc.SetParameter(0,1)
   #      fitFunc.SetParLimits(0,3,2)
   #      
   #   fitFunc.SetParameter(1,0)
   #   fitFunc.SetParLimits(1,3,2)
      
      
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


