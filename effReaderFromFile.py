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
#from helpers import *

def drange(start, stop, step):
  r = start
  while r < stop:
    yield r
    r += step


#######################################

class effReaderFromFile:
   
   def __init__(self,folder,process,benergy,cat,massLow=115.,massHigh=160.):
      # this is the folder where all the extrapolation
      # parameters are contained
      self.folder   = folder
      # process should be GluGlu or VBF
      self.process  = process

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
      self.fitFuncs ['eff'] = root.TF1("eff_%s"  % cat, "pol2",massLow,massHigh)

      self.massLow  = massLow
      self.massHigh = massHigh


   def initPol2(self,function,process,benergy,category,par):

      filename = self.folder + "/effExtrapolation_" + process + '_' + benergy + '_' + category + '.txt'

      file = open(filename)
      for line in file:
         if ('#' in line):
            continue
         else:
            line = line.split('\n')[0]
            line = line.split(' ')
            #print line

            # the function is a pol2 
            for id in range(0,3):
               parValue = line[(2*id)]
               parError = line[(2*id)+1]
               
               function.SetParameter( int(id),float(parValue) )
               function.SetParError ( int(id),float(parError) )
                                     
      
   def initFunctions(self):

      # for the extrapolations we used pol2 functions and the parameters
      self.initPol2(self.fitFuncs['eff'], self.process,self.benergy,self.category, 'eff' )


   def getEff(self):

      self.initPol2(self.fitFuncs['eff'], self.process,self.benergy,self.category, 'eff' )
      func = self.fitFuncs['eff']
      
      efficiency = {}

      massrange = drange(self.massLow, self.massHigh+0.5, 0.5)
      for mass in massrange:
         parValue = func.Eval(mass)
         #print parname, mass, parValue
         efficiency['%s' % mass] = float(parValue)

      return efficiency

     
      
#######################################

#if __name__ == "__main__":

#  effReader = effReaderFromFile('fitresults',
#                                'gg',
#                                '8TeV',
#                                'Jets01PassPtG10BO')
#  
#  # example with one parameter
#  efficiencies = effReader.getEff()
#  massrange = drange(115.0, 150.5, 0.5)
#  for mass in massrange:
#    print mass, efficiencies['%s' % mass]

