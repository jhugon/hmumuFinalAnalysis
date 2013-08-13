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



#######################################

class signalPars:
   
   def __init__(self,folder,process,benergy,cat,massLow=115.,massHigh=160.,resSF=1.):
      # this is the folder where all the extrapolation
      # parameters are contained
      self.folder   = folder
      # process should be GluGlu or VBF
      self.process  = process
      self.resSF = resSF

      # this is to match Justin's definitions without
      # obliging him to rewrite his code
      if (process == 'GluGlu'):
         self.process = 'gg'
      if (process == 'VBF'):
         self.process = 'vbf'

      # benergy should be 7TeV or 8TeV
      if benergy == "14TeV":
        benergy = "8TeV"
      self.benergy  = benergy
      # the current categories available are:
      #    - GluGlu BB,BO,BE,OO,OE,EE
      #    - VBF, VBFCiCLoose, VBFCiCTight
      self.category = cat

      self.fitFuncs = {}
      self.fitFuncs ['meanG1' ] = root.TF1("meanG1_%s"  % cat, "pol1",massLow,massHigh)
      self.fitFuncs ['widthG1'] = root.TF1("widthG1_%s" % cat, "pol1",massLow,massHigh)
      self.fitFuncs ['meanG2' ] = root.TF1("meanG2_%s"  % cat, "pol1",massLow,massHigh)
      self.fitFuncs ['widthG2'] = root.TF1("widthG2_%s" % cat, "pol1",massLow,massHigh)
      self.fitFuncs ['mixGG'  ] = root.TF1("mixGG_%s"   % cat, "pol1",massLow,massHigh)

      self.massLow  = massLow
      self.massHigh = massHigh

      # useful members to draw the signal templates
      self.canvas = None
      self.sigTemplates = {}

      # title and savename for the signal templates
      self.title = "CMS Simulation "
      self.savename  = folder + '/png/sigTemplates_'
      self.savename += process + '_'
      self.savename += benergy + '_'
      self.savename += cat + '_DG.png'

      if (process == 'gg' or process == 'GluGlu'):
         self.title += 'gg #rightarrow H #rightarrow #mu#mu at 8 TeV'

         if self.category == 'BB':
            self.title += ', Barrel-Barrel'
         if self.category == 'BO':
            self.title += ', Barrel-Overlap'
         if self.category == 'BE':
            self.title += ', Barrel-Endcap'
         if self.category == 'OO':
            self.title += ', Overlap-Overlap'
         if self.category == 'OE':
            self.title += ', Overlap-Endcap'
         if self.category == 'EE':
            self.title += ', Endcap-Endcap'

      if (process == 'vbf' or process == 'VBF'):
          self.title += 'VBF H #rightarrow #mu#mu at 8 TeV'
          
          if self.category == 'VBF':
             self.title += ', BDT'
          if self.category == 'VBFCiCLoose':
             self.title += ', Loose'
          if self.category == 'VBFCiCTight':
             self.title += ', Tight'

      
   def initPol1(self,function,process,benergy,category,par):
      if benergy == "14TeV":
        benergy = "8TeV"

      filename = self.folder + "/extrapolation_" + par + '_' + process + '_' + benergy + '_' + category + '_DG.txt'

      file = open(filename)
      for line in file:
         if ('#' in line):
            continue
         else:
            line = line.split('\n')[0]
            line = line.split(' ')
            #print line

            # the function is a pol1 
            for id in range(0,2):
               parValue = line[(2*id)]
               parError = line[(2*id)+1]
               
               function.SetParameter( int(id),float(parValue) )
               function.SetParError ( int(id),float(parError) )
                                     
      
   def initFunctions(self):

      # for the extrapolations we used pol1 functions and the parameters
      # of the double Gaussian are 5
      self.initPol1(self.fitFuncs['meanG1' ], self.process,self.benergy,self.category, 'meanG1' )
      self.initPol1(self.fitFuncs['widthG1'], self.process,self.benergy,self.category, 'widthG1')
      self.initPol1(self.fitFuncs['meanG2' ], self.process,self.benergy,self.category, 'meanG2' )
      self.initPol1(self.fitFuncs['widthG2'], self.process,self.benergy,self.category, 'widthG2')
      self.initPol1(self.fitFuncs['mixGG'  ], self.process,self.benergy,self.category, 'mixGG'  )


   def getPar(self,parname):

      self.initPol1(self.fitFuncs[parname], self.process,self.benergy,self.category, parname )
      func = self.fitFuncs[parname]
      
      parameter = {}

      massrange = drange(self.massLow, self.massHigh+0.5, 0.5)
      for mass in massrange:
         parValue = func.Eval(mass)
         #print parname, mass, parValue
         parameter['%s' % mass] = float(parValue)
         if parname == "widthG1":
           parameter['%s' % mass] = parameter['%s' % mass]*self.resSF

      return parameter


   # this method returns all the parameters
   def getPars(self):

      par_meanG1  = self.getPar( 'meanG1')
      par_widthG1 = self.getPar('widthG1')
      par_meanG2  = self.getPar( 'meanG2')
      par_widthG2 = self.getPar('widthG2')
      par_mixGG   = self.getPar(  'mixGG')

      return par_meanG1, par_widthG1, par_meanG2, par_widthG2, par_mixGG


   # draw the points every X GeV step
   def draw(self, step):

      minMass   = self.massLow -10.0
      maxMass   = self.massHigh+10.0
      mMuMu     = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)


      # get the parameters
      par_meanG1, par_widthG1, par_meanG2, par_widthG2, par_mixGG = self.getPars()
      
      massrange = drange(self.massLow, self.massHigh+0.5, step)
      #workspace = root.RooWorkspace("w")
      #wImport = getattr(workspace,"import")
      self.canvas = root.TCanvas("canvas","",0,0,1100,650)
      self.canvas.cd()
      plotMmumu = mMuMu.frame(minMass,maxMass)
      
      for mass in massrange:

         #if (mass < 140):
         #   continue
      
         #####################################################################
         # define the Double Gaussian
         meanG1 = root.RooRealVar("MeanG1_%s"%mass,"MeanG1_%s"%mass, par_meanG1['%s'%mass], 100.,170.)
         meanG2 = root.RooRealVar("MeanG2_%s"%mass,"MeanG2_%s"%mass, par_meanG2['%s'%mass], 100.,170.)
 
         widthG1 = root.RooRealVar("WidthG1_%s"%mass,"WidthG1_%s"%mass, par_widthG1['%s'%mass], 0.1,20.0)
         widthG2 = root.RooRealVar("WidthG2_%s"%mass,"WidthG2_%s"%mass, par_widthG2['%s'%mass], 0.1,5.0)
 
         mixGG = root.RooRealVar("mixGG_%s"%mass,"mixGG_%s"%mass, 1-par_mixGG['%s'%mass], 0.0,1.0)
         
         gaus1 = root.RooGaussian("gaus1_%s"%mass,"gaus1_%s"%mass,mMuMu,meanG1,widthG1)
         #gaus2 = root.RooGaussian("gaus2_%s"%mass,"gaus2_%s"%mass,mMuMu,meanG2,widthG2)
         gaus2 = root.RooGaussian("gaus2_%s"%mass,"gaus2_%s"%mass,mMuMu,meanG1,widthG2)
 
         pdfMmumuGG = root.RooAddPdf("pdfMmumuGG_%s"%mass,"pdfMmumuGG_%s"%mass,gaus1,gaus2,mixGG)

         #wImport(pdfMmumuGG)
         
         print mass
         self.sigTemplates['%s' % mass] = pdfMmumuGG

         plotMmumu.GetXaxis().SetTitle("M(#mu#mu) [GeV/c^{2}]")
         plotMmumu.GetYaxis().SetTitle("a.u.")
         plotMmumu.SetTitle(self.title)
         pdfMmumuGG.plotOn(plotMmumu,root.RooFit.LineColor(root.kAzure+2))
         plotMmumu.Draw()

         self.canvas.Update()

      self.canvas.SaveAs(self.savename)
      
      
#######################################

if __name__ == "__main__":

   s = signalPars('fitresults',
                  'gg',
                  '8TeV',
                  'Jets01PassPtG10BB')
                  #'Jets01PassPtG10BO')
                  #'Jets01PassPtG10BE')
                  #'Jets01FailPtG10BB')
                  #'Jets01FailPtG10BO')
                  #'Jets01FailPtG10BE')
                  #'Jets01PassPtG10EE')
                  #'Jets01FailPtG10EE')

   #example with one parameter
   #parameter = s.getPar('meanG1')
   #massrange = drange(115.0, 150.5, 0.5)
   #for mass in massrange:
   #   print mass, parameter[ '%s' % mass ]
      
   # example with all the parameters
   #par_meanG1, par_widthG1, par_meanG2, par_widthG2, par_mixGG = s.getPars()
   #for mass in par_meanG1.keys():
   #   print mass, par_meanG1[mass], par_widthG1[mass], par_meanG2[mass], par_widthG2[mass], par_mixGG[mass] 

   # example to draw signal templates every 5 GeV/c2
   s.draw(5)


#
#   s = signalPars('fitresults',
#                  'gg',
#                  '8TeV',
#                  'IncPreselPtG10BB')
#
#   #example with one parameter
#   parameter = s.getPar('meanG1')
#   massrange = drange(115.0, 150.5, 0.5)
#   for mass in massrange:
#      print mass, parameter[ '%s' % mass ]
#      
#   # example with all the parameters
#   #par_meanG1, par_widthG1, par_meanG2, par_widthG2, par_mixGG = s.getPars()
#   #for mass in par_meanG1.keys():
#   #   print mass, par_meanG1[mass], par_widthG1[mass], par_meanG2[mass], par_widthG2[mass], par_mixGG[mass] 
#
#   # example to draw signal templates every 5 GeV/c2
#   s.draw(5)

