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

#folder = "/afs/cern.ch/user/d/digiovan/public/forJustin/signalfits7and8TeV/fitresults/"
folder = "signalfits/fitresults/"

baseNames = ["IncPreselPtG10BB",
             "IncPreselPtG10BO",
             "IncPreselPtG10BE",
             "IncPreselPtG10OO",
             "IncPreselPtG10OE",
             "IncPreselPtG10EE",
             "VBFBDTCut"]

fitNames= ['DG','CBG']
benergy = ['7TeV','8TeV']


##########################################################################################
def expandParameters(baseName,fit,benergy):

   parameters = {}
   parameters[benergy] ={}

   listOfFiles = commands.getoutput('ls %s | grep %s | grep %s | grep %s | grep fit'
                                    % (folder,baseName,benergy,fit) ).split('\n')
   #print listOfFiles
   
   masses = []
   for file in listOfFiles:
      input = open(folder+"/"+file)
   
      for line in input:
         if ('#' in line):
            continue
         else:
            line = line.split('\n')[0]
            line = line.split(' ')
            #print line
   
            mass = line[1]
            masses.append( mass)
   
   for mass in masses:
      # add the mass
      parameters[benergy][mass] = {}
      
      # add categories, IncPreselPtG10BB, IncPreselPtG10BO, etc, etc
      parameters[benergy][mass][baseName] = {}

      varNames    = ['meanG1','widthG1','meanG2','widthG2','mixGG']
      varNamesEE  = ['meanSG','widthSG']
      varNamesCBG = ['mean','width1','width2','mix','Alpha','n']
   
      if ("EE" in baseName):
         varNames = varNamesEE
   
      if (fit == 'CBG'):
         varNames = varNamesCBG

      #print varNames

   
      # add the fit name
      parameters[benergy][mass][baseName][fit] = {}
   
   
      # add the fit results
      for file in listOfFiles:
         if ( mass not in file ):
            continue
         #print file
   
         # open the file with the par values and errors
         input = open(folder+"/"+file)
         
         for line in input:
            if ('#' in line):
               continue
            else:
               line = line.split('\n')[0]
               line = line.split(' ')
               #print line
    
            for id in range(0,len(varNames)):
               parValue = line[2+(2*id)]
               parError = line[2+(2*id)+1]
   
               parameters[benergy][mass][baseName][fit][varNames[id]]          = float(parValue)
               parameters[benergy][mass][baseName][fit][varNames[id]+'_error'] = float(parError)

   return parameters
   
##########################################################################################
   
# this is the dictionary we want to fill
parameters = {}



for baseName in baseNames:
   print  baseName

   for fit in fitNames:
      #print fit

      for bener in benergy:
         parTmp = expandParameters(baseName,fit,bener)

         if (len(parameters.keys()) == 0 ):
             parameters = parTmp
             
         for ben in parTmp.keys():
            #print ben, parameters.has_key(ben)
            if ( parameters.has_key(ben) == False ):
               parameters[ben] = parTmp[ben]
               
         
            for mass in parTmp[ben].keys():
               #print mass, parameters[ben].has_key(mass)
         
               if ( parameters[ben].has_key(mass) == False ):
                  parameters[ben][mass] = parTmp[ben][mass]
         
         
               for cat in parTmp[ben][mass].keys():
                  #print cat, parameters[ben][mass].has_key(cat)
         
                  if ( parameters[ben][mass].has_key(cat) == False ):
                     parameters[ben][mass][cat] = parTmp[ben][mass][cat]
         
         
                  for thefit in parTmp[ben][mass][cat].keys():
                     #print thefit
         
                     if ( parameters[ben][mass][cat].has_key(thefit) == False ):
                        parameters[ben][mass][cat][thefit] = parTmp[ben][mass][cat][thefit]
                        

