#!/usr/bin/env python

import sys
import os
import commands
from helpers import * 
import time

# ------------------------------------------------------------------------------

def Usage():
   print 'Wrong syntax: \n'
   print '   ./resubmitMassLimitlxbatch [dir (/afs/cern.ch/work/d/digiovan/code/CMSSW_6_1_1/batchBaselinePP/)]\n'
   sys.exit()
    
# ------------------------------------------------------------------------------

# ====================   
# set useful variables
if (len(sys.argv) > 2 ):
   Usage()

# default values
outDir = '/afs/cern.ch/work/d/digiovan/code/CMSSW_6_1_1/batchBaselinePP/'
if (len(sys.argv) == 2):
   outDir  =  sys.argv[1]

print '\n----------------------------------------------------------------------'
print 'ANALYZING OUTPUTS IN DIR %s' % outDir
# ====================

PWD = os.getenv('PWD')

total = 0
done  = 0
toResubmit = 0

resubScriptName = '%s/resubmit.sh' % outDir 
resubScript = open(resubScriptName,'w') 
resubScript.write('#! /bin/bash\n\n')

massrange = drange(115, 160.5, 1)
for mass in massrange:

  outfiles   = commands.getoutput( "ls %s | grep out | grep %s"  % (outDir,mass) ).split("\n") 
  statsCards = commands.getoutput( "ls statsCards | grep txt | grep %s" % (mass) ).split("\n") 

  # remove the ".out" ending to outfiles
  outfiles = [file.replace(".out","") for file in outfiles]

  # compute the difference
  diffList=list(set(statsCards)-set(outfiles))

  #print diffList

  total      += int ( len(statsCards) )
  done       += int ( len(outfiles)   )
  toResubmit += int ( len(diffList)   )

  # write the jobs to resubmit
  for card in diffList:
    resubScript.write("bsub -q cmscaf1nh -o /dev/null lxbatch.sh %s\n" % card )


resubScript.close()

print "total # of events: %s" % total
print "done # of events: %s" %  done
print "resubmitted # of events: %s" % toResubmit

# this gives some time to realize how many jobs are going to be submitted
print "going to resubmit them... wait"
time.sleep(5)

resubCommands  = " cd %s; " % outDir
resubCommands += "eval \`scramv1 runtime -sh\`; "
resubCommands += "bash resubmit.sh; "
resubCommands += "bash getStatus2.sh; "
resubCommands += "cd %s" % PWD

copyCommands  = "echo \"Copying output files from lxplus5...\"; "
copyCommands += "scp lxplus5:%s/*.out statsInput/.; " % outDir
copyCommands += "scp lxplus5:%s/*.sig statsInput/.; " % outDir
copyCommands += "scp lxplus5:%s/*.mu statsInput/.; " % outDir
copyCommands += "scp lxplus5:%s/*.txt.root statsInput/.; " % outDir
copyCommands += "scp lxplus5:%s/*.png statsInput/.; " % outDir
copyCommands += "scp lxplus5:%s/*CCC*.root statsInput/.; " % outDir

print resubCommands
os.system( resubCommands )

print copyCommands
os.system( copyCommands )

print '----------------------------------------------------------------------\n'

