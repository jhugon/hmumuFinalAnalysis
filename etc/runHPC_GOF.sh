#!/bin/bash

NJOBS=`ls *.txt | wc -l`
sed "s/YAYYAYYAY/$NJOBS/" gofHPC_Template.sh > gof.sh

chmod +x gof.sh

qsub gof.sh 1>&2

echo "Submitted $NJOBS jobs" 1>&2

echo "now working on the background only jobs..." 1>&2

read -d '' pythonBakListCommand <<"EOF"
import glob
import re
files = glob.glob("*.txt")
cats = set()
for f in files:
  m = re.match(r"(.+_.+TeV)_([0-9.]+)\.txt",f)
  cat =  m.group(1)
  if not cat in cats:
    cats.add(cat)
result = 0
for cat in cats:
  files = glob.glob(cat+"*.txt")
  for f in files:
    result += 1
    break
print result
EOF
BACKGROUNDN=`python -c "$pythonBakListCommand"`

sed "s/YAYYAYYAY/$BACKGROUNDN/" gofHPC_Template_bak.sh > gof_bak.sh

chmod +x gof_bak.sh

qsub gof_bak.sh 1>&2

echo "Submitted $BACKGROUNDN Bak jobs" 1>&2
echo $BACKGROUNDN
