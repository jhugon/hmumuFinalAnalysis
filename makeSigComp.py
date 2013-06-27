#!/usr/bin/env python

from helpers import *
import ROOT as root
import glob
import re
import os.path
from copy import deepcopy
from array import array

from xsec import *
import commands

def setHist(histo, color):
   histo.SetFillStyle(3001)
   histo.SetFillColor(color)
   histo.SetLineColor(color)

# ------------------------------------------------------------------------------

def Usage():
   print 'Wrong syntax: \n'
   print '   ./makeSigComp.py [mass] [period (7TeV,8TeV)]\n'
   sys.exit()
    
# ------------------------------------------------------------------------------

# ====================   
# set useful variables
if (len(sys.argv) < 3):
   Usage()

# default values
mass   = sys.argv[1]
period = sys.argv[2]


setStyle()

dirName = "statsCards/"

allfiles = ['Jets01PassCatAll',
            'Jets01FailCatAll',
            'Jet2CutsVBFPass',
            'Jet2CutsGFPass',
            'Jet2CutsFailVBFGF'
            ]

#~48 Charactars Max
titleMap = {

  "Jets01PassPtG10BB": "Non-VBF Tight BB",
  "Jets01PassPtG10BO": "Non-VBF Tight BO",
  "Jets01PassPtG10BE": "Non-VBF Tight BE",
  "Jets01PassPtG10OO": "Non-VBF Tight OO",
  "Jets01PassPtG10OE": "Non-VBF Tight OE",
  "Jets01PassPtG10EE": "Non-VBF Tight EE",
  "Jets01PassCatAll" : "Non-VBF Tight",
                               
  "Jets01FailPtG10BB": "Non-VBF Loose BB",
  "Jets01FailPtG10BO": "Non-VBF Loose BO",
  "Jets01FailPtG10BE": "Non-VBF Loose BE",
  "Jets01FailPtG10OO": "Non-VBF Loose OO",
  "Jets01FailPtG10OE": "Non-VBF Loose OE",
  "Jets01FailPtG10EE": "Non-VBF Loose EE",
  "Jets01FailCatAll" : "Non-VBF Loose",
                               
  "Jets01SplitCatAll": "Non-VBF",


  "Jet2CutsVBFPass"  :"VBF, VBF Tight",
  "Jet2CutsGFPass"   :"VBF, GF Tight",
  "Jet2CutsFailVBFGF":"VBF, Loose",

  "Jet2SplitCutsGFSplit" : "VBF Preselection",
  "CombSplitAll" : "Combination",

}


# done in order not to have all the bins attached to each other
nbins = 2*len(allfiles)+5

# Absolute yield
hSig_ggH = root.TH1F("hSig_ggH", "", nbins, 0, nbins)
hSig_qqH = root.TH1F("hSig_qqH", "", nbins, 0, nbins)

# Percentage yield
hPerc_ggH = root.TH1F("hPerc_ggH", "", nbins, 0, nbins)
hPerc_qqH = root.TH1F("hPerc_qqH", "", nbins, 0, nbins)



# var to fill the correct bin (every 2 bins, just for
# graphical needs) 
iBin = 2

# loop
for fn in allfiles:

   # get the card
   card = dirName + fn + '_' + period + '_' + mass + '.txt'
   #print card

   # get the # of events
   line = commands.getoutput('grep rate %s' % card)

   rates = line.split()
   #print rates
   
   # fill the counters
   ggH = 0.
   qqH = 0.
   
   nEntries = int( (len(rates)-1)/3. )
   for i in range(1,nEntries+1):

      id_ggH = (i - 1)*3 + 1
      id_qqH = (i - 1)*3 + 2
      
      ggH += float(rates[id_ggH])
      qqH += float(rates[id_qqH])

      #print  id_ggH, id_qqH, float(rates[id_ggH]), float(rates[id_qqH])

   signal = ggH + qqH
   #print iBin, ggH, qqH, signal, ggH/signal, qqH/signal
   
   hSig_ggH.SetBinContent(iBin,ggH)
   hSig_qqH.SetBinContent(iBin,qqH)    

   hPerc_ggH.SetBinContent(iBin,ggH/signal)
   hPerc_qqH.SetBinContent(iBin,qqH/signal)    

   iBin += 2



# set some nice style for plots
setHist(hSig_ggH,  629)
setHist(hPerc_ggH, 629)

setHist(hSig_qqH,  root.kAzure+2)
setHist(hPerc_qqH, root.kAzure+2)

# create the stacks
stackSig = root.THStack("stackSig","CMS Simulation, H #rightarrow #mu #mu, m_{H}=%s GeV/c^{2} at %s" % (mass, period))
stackSig.Add( hSig_ggH )
stackSig.Add( hSig_qqH )

stackPerc = root.THStack("stackPerc","CMS Simulation, H #rightarrow #mu #mu, m_{H}=%s GeV/c^{2} at %s" % (mass, period))
stackPerc.Add( hPerc_ggH )
stackPerc.Add( hPerc_qqH )

# TLegend
#urLegendPos = [0.65,0.67,0.9,0.9]
urLegendPos = [0.76,0.70,0.93,0.85]
leg = root.TLegend(*urLegendPos)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.042)
leg.AddEntry(hSig_ggH," GF", "f")
leg.AddEntry(hSig_qqH," VBF", "f")
# draw the stacks

# absolute yield
canvasSig = root.TCanvas("canvasSig","",0,0,800,600)
canvasSig.cd()

# I need to draw it first otherwise its pointer
# will not exist and setting the binlabel will not work...
stackSig.Draw()

stackSig.GetYaxis().SetTitle("Expected Yield")
stackSig.GetXaxis().SetTitle("Category")

#print nbins
for id in range(0, len(allfiles)):

   bin = 2*id + 2
   #print bin
   stackSig.GetXaxis().SetBinLabel( bin, titleMap[allfiles[id]] )

stackSig.Draw()

leg.Draw("same")
canvasSig.Update()
canvasSig.SaveAs( "plots/png/yieldAbs_m%s_%s.png"  % (mass,period) )
canvasSig.SaveAs("plots/root/yieldAbs_m%s_%s.root" % (mass,period) )


# relative yield
canvasPerc = root.TCanvas("canvasPerc","",900,0,800,600)
canvasPerc.cd()

# I need to draw it first otherwise its pointer
# will not exist and setting the binlabel will not work...
stackPerc.Draw()

stackPerc.GetYaxis().SetTitle("Excted Percentage Yield")
stackPerc.GetXaxis().SetTitle("Category")

for id in range(0, len(allfiles)):
   bin = 2*id + 2
   stackPerc.GetXaxis().SetBinLabel( bin, titleMap[allfiles[id]] )

stackPerc.Draw()

leg.Draw("same")
canvasPerc.Update()
canvasPerc.SaveAs( "plots/png/yieldPerc_m%s_%s.png"  % (mass,period) )
canvasPerc.SaveAs("plots/root/yieldPerc_m%s_%s.root" % (mass,period) )

