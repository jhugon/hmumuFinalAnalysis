#! /usr/bin/env python

import sys
import os
import commands
from helpers import * 

# ------------------------------------------------------------------------------

def Usage():
   print 'Wrong syntax: \n'
   print '   ./makeAllCardsLxbatch.py [queue_type(cmscaf1nd)]\n'
   sys.exit()
    
# ------------------------------------------------------------------------------

# ====================   
# set useful variables
if (len(sys.argv) > 2 ):
   Usage()

# default values
queue = 'cmscaf1nd'
if (len(sys.argv) == 2):
   queue  =  sys.argv[1]

print 'SUBMIT ON BATCH QUEUE %s' % queue
print '----------------------------------------------------------------------\n'
# ====================

PWD = os.getenv('PWD')
massrange = drange(115, 155.5, 1)
#massrange = drange(125, 127.5, 1)
for mass in massrange:

   
   jobscript = 'JOBS/' + 'jobs_makeCards_mass%s.sh' % str(mass).replace(".","p")
   outfile = open(jobscript,'w') 
   outfile.write('#! /bin/bash\n\n')
   outfile.write('cd /afs/cern.ch/work/d/digiovan/code/CMSSW_6_1_1/src/\n')
   outfile.write('source ../STARTUP\n')
   outfile.write('cd %s\n' % PWD)
   outfile.write('./makeCards.py -m %s\n\n' % mass)
   outfile.close()
   
   print '\nbsub -q %s -J JOB_%s < %s' \
         % (queue, mass, jobscript)
   
   os.system( '\nbsub -q %s -J JOB_%s < %s' \
              % (queue, mass, jobscript)
              )
print '----------------------------------------------------------------------\n\n'


