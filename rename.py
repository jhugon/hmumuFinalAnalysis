#!/usr/bin/env python

import sys
import os
import commands


# ------------------------------------------------------------------------------

def Usage():
   print 'Wrong syntax: \n'
   print '   ./rename.py GluGlu 8TeV\n'
   sys.exit()
    
# ------------------------------------------------------------------------------

# ====================   
# set useful variables
if (len(sys.argv) > 3 ):
   Usage()

# default values
process = 'GluGlu'
benergy = '8TeV'

if (len(sys.argv) == 2):
   process  =  sys.argv[1]

if (len(sys.argv) == 3):
   process  =  sys.argv[1]
   benergy  =  sys.argv[2]
   
print '\n----------------------------------------------------------------------'
print 'RENAMING FILES WITH PROCESS %s AT %s' % (process, benergy)
# ====================

PWD = os.getenv('PWD')

#Eff8TeV_ggHiggs125.txt Eff8TeV_vbfHiggs125.txt
endprocess = 'gg'
if (process == "VBF"):
   endprocess = 'vbf'
   
list = commands.getoutput( "ls | grep %s | grep %s"  % (process,benergy) ).split("\n") 

for file in list:

   mass = file.split("_")[0].replace("%sHmumu"%process,"").replace("p5",".5")
   print "cp %s Eff%s_%sHiggs%s.txt" % (file,benergy,endprocess,mass)
   os.system("cp %s Eff%s_%sHiggs%s.txt" % (file,benergy,endprocess,mass))

