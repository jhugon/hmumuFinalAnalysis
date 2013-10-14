#!/usr/bin/env python

from scipy import *
import matplotlib.pyplot as mpl
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

from singleHelpers import *
from signalPars import *

TITLEMAP = {
  "Jets01PassPtG10BB": "0,1-Jet Tight BB",
  "Jets01PassPtG10BO": "0,1-Jet Tight BO",
  "Jets01PassPtG10BE": "0,1-Jet Tight BE",
  "Jets01PassPtG10OO": "0,1-Jet Tight OO",
  "Jets01PassPtG10OE": "0,1-Jet Tight OE",
  "Jets01PassPtG10EE": "0,1-Jet Tight EE",
  #"Jets01PassCatAll" : "0,1-Jet Tight Combination",
                        
#  "Jets01FailPtG10BB": "0,1-Jet Loose BB",
#  "Jets01FailPtG10BO": "0,1-Jet Loose BO",
#  "Jets01FailPtG10BE": "0,1-Jet Loose BE",
#  "Jets01FailPtG10OO": "0,1-Jet Loose OO",
#  "Jets01FailPtG10OE": "0,1-Jet Loose OE",
#  "Jets01FailPtG10EE": "0,1-Jet Loose EE",
#  #"Jets01FailCatAll" : "0,1-Jet Loose Combination",

  "Jet2CutsVBFPass":"2-Jet VBF Tight",
  "Jet2CutsGFPass":"2-Jet GF Tight",
  "Jet2CutsFailVBFGF":"2-Jet Loose",
}

categories = sorted(TITLEMAP.keys())

energies = ["7TeV","8TeV"]

hmasses = range(115,156)

def annotationsTop(ax,annotation):
  ax.text(0,1.03, 'CMS Preliminary',
        size = 'large',
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax.transAxes)
  ax.text(1,1.03, annotation,
        size = 'large',
        horizontalalignment='right',
        verticalalignment='bottom',
        transform=ax.transAxes)

fig = mpl.figure()
fig2 = mpl.figure()
for energy in energies:
  ax1 = fig.add_subplot(111)
  ax1.set_xlabel(r"$m_H$ $[GeV/c^2]$")
  ax1.set_ylabel(r"Narrow Gaussian $\sigma$ $[GeV/c^2]$")
  ax1.set_xlim(110,160)
  ax1.set_ylim(1.,5.)
  ax2 = fig2.add_subplot(111)
  ax2.set_xlabel(r"$m_H$ $[GeV/c^2]$")
  ax2.set_ylabel(r"$5\times$ Narrow Gaussian Width $[GeV/c^2]$")
  ax2.set_xlim(110,160)
  ax2.set_ylim(10.,50.)
  plots = []
  for cat in categories:
    catTitle = TITLEMAP[cat]
    prodMode = 'gg'
    linestyle = '-'
    if 'Jet2' in cat:
      prodMode = 'vbf'
      linestyle = '--'
    sigFits = signalPars('fitresults',prodMode,energy,cat)
    sigmaDict = sigFits.getPar('widthG1')
    sigmas = [sigmaDict["{0:.1f}".format(hmass)] for hmass in hmasses]
    sigmas5 = [10*x for x in sigmas]
    tmpPlot = ax1.plot(hmasses,sigmas,ls=linestyle,label=catTitle)
    tmpPlot2 = ax2.plot(hmasses,sigmas5,ls=linestyle,label=catTitle)
  handles, labels = ax1.get_legend_handles_labels()
  handles2, labels2 = ax2.get_legend_handles_labels()
  ax1.legend(handles, labels,prop={'size':'small'},ncol=2,loc='upper left',frameon=False)
  ax2.legend(handles2, labels2,prop={'size':'small'},ncol=2,loc='upper left',frameon=False)
  annotationsTop(ax1,'$\sqrt{s}='+energy[0]+'$ TeV')
  annotationsTop(ax2,'$\sqrt{s}='+energy[0]+'$ TeV')
  fig.savefig("SignalWidthVMH_"+energy+".png")
  fig2.savefig("Signal5WidthVMH_"+energy+".png")
  fig.clf()
  fig2.clf()
