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
   logfile = 'JOBS/' + 'log_makeCards_mass%s.sh' % str(mass).replace(".","p")
   outfile = open(jobscript,'w') 
   outfile.write('#! /bin/bash\n\n')
   outfile.write("""

echo "Sourcing cmsset_default.sh"
cd /afs/cern.ch/cms
source cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc472
echo "SCRAM_ARCH is $SCRAM_ARCH"
cd $LS_SUBCWD
echo "In Directory: "
pwd
eval `scramv1 runtime -sh`
echo "cmsenv success!"
date

"""
)
   outfile.write('./makeCards.py -m %s\n\n' % mass)
   outfile.write('date\n\n')
   outfile.close()
   
   print '\nbsub -q %s -J JOB_%s -o %s < %s' \
         % (queue, mass, logfile, jobscript)
   
   os.system( '\nbsub -q %s -J JOB_%s -o %s < %s' \
              % (queue, mass, logfile, jobscript)
              )
print '----------------------------------------------------------------------\n\n'


