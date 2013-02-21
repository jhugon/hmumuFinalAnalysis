#!/usr/bin/env python

import sys
import ROOT

ROOT.gROOT.SetBatch(True)

if sys.argv < 2:
  print("Error: Requires filename argument")
  sys.exit(1)

file = ROOT.TFile(sys.argv[1])
if file == None: raise RuntimeError, "Cannot open file %s" % "mlfit.root"
fit_s  = file.Get("fit_s")
fit_b  = file.Get("fit_b")
prefit = file.Get("nuisances_prefit")
if fit_s == None or fit_s.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the signal fit 'fit_s'"     % args[0]
if fit_b == None or fit_b.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the background fit 'fit_b'" % args[0]
if prefit == None or prefit.ClassName() != "RooArgSet":    raise RuntimeError, "File %s does not contain the prefit nuisances 'nuisances_prefit'"  % args[0]

table = {}
fpf_b = fit_b.floatParsFinal()
fpf_s = fit_s.floatParsFinal()
for i in range(fpf_s.getSize()):
    nuis_s = fpf_s.at(i)
    name   = nuis_s.GetName();
    nuis_b = fpf_b.find(name)
    nuis_p = prefit.find(name)
    row = []
    flag = False;
    # .getVal() .getError() .getMin() .getMax()

    rowString = ""
    rowString += "{0:30}".format(name)
    rowString += "{0:8.3f} +/- {1:5.3f}  ".format(nuis_s.getVal(),nuis_s.getError())

    rowString += "{0:15}".format("[{0:.1f},{1:.1f}]".format(nuis_s.getMin(),nuis_s.getMax()))

    if nuis_p:
      rowString += " Gaus Constraint: {0:8.3f} +/- {1:5.3f}".format(nuis_p.getVal(),nuis_p.getError())

    
    
    print(rowString)
