#!/usr/bin/env python

from xsec import *
from helpers import *
import ROOT as root
import os
import sys
from math import ceil

mass = '125'

energies=[
'7TeV',
'8TeV'
]

massCut = "dimuonMass > 110. && dimuonMass < 170."
jet01Cuts = " && !(jetLead_pt>40. && jetSub_pt>30. && ptMiss<40.)"
jet2Cuts = " && jetLead_pt>40. && jetSub_pt>30. && ptMiss<40."
errors = ["JES","JER","PUID","MCStat"]

effReader = EfficiencyReader()

datasets = ["gg","vbf"]
categories = [
  "Jets01PassPtG10",
  "Jets01FailPtG10",
  "Jet2CutsVBFPass",
  "Jet2CutsGFPass",
  "Jet2CutsFailVBFGF"
]
labels = [
  "Non-VBF Presel. Tight",
  "Non-VBF Presel. Loose",
  "VBF-Presel. VBF Tight",
  "VBF-Presel. GF Tight",
  "VBF-Presel. Loose"
]

class ErrorStruct(object):
  def __init__(self):
    self.contents = []
  def __setattr__(self,name,value):
    if hasattr(self,"contents"):
      #print("In my overloaded function!! name: %s value: %s" %(name,value))
      #if len(self.contents) < 5:
        self.contents.append(name)
    super(ErrorStruct,self).__setattr__(name,value)
  def __getitem__(self,key):
    return self.__getattribute__(key)

errors = ErrorStruct()

######################################
######################################
######################################
######################################
######################################

errors.MCStat = {
  'gg' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : 1.0244,
      'Jet2CutsGFPass' : 1.0485,
      'Jet2CutsVBFPass' : 1.1019,
      'Jets01FailPtG10' : 1.0110,
      'Jets01PassPtG10' : 1.0042,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : 1.0222,
      'Jet2CutsGFPass' : 1.0429,
      'Jet2CutsVBFPass' : 1.0732,
      'Jets01FailPtG10' : 1.0109,
      'Jets01PassPtG10' : 1.0043,
      },
    },
  'vbf' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : 1.0136,
      'Jet2CutsGFPass' : 1.0115,
      'Jet2CutsVBFPass' : 1.0110,
      'Jets01FailPtG10' : 1.0447,
      'Jets01PassPtG10' : 1.0052,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : 1.0127,
      'Jet2CutsGFPass' : 1.0112,
      'Jet2CutsVBFPass' : 1.0096,
      'Jets01FailPtG10' : 1.0456,
      'Jets01PassPtG10' : 1.0054,
      },
    },
  'w' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    },
  'z' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    },
  }

errors.JES = {
  'gg' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : 1.0832,
      'Jet2CutsGFPass' : 1.0585,
      'Jet2CutsVBFPass' : 1.0759,
      'Jets01FailPtG10' : -1.0006,
      'Jets01PassPtG10' : -1.0054,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : 1.0799,
      'Jet2CutsGFPass' : 1.0477,
      'Jet2CutsVBFPass' : 1.0606,
      'Jets01FailPtG10' : -1.0006,
      'Jets01PassPtG10' : -1.0060,
      },
    },
  'vbf' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : 1.0572,
      'Jet2CutsGFPass' : 1.0381,
      'Jet2CutsVBFPass' : 1.0358,
      'Jets01FailPtG10' : -1.0321,
      'Jets01PassPtG10' : -1.0308,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : 1.0482,
      'Jet2CutsGFPass' : 1.0249,
      'Jet2CutsVBFPass' : 1.0185,
      'Jets01FailPtG10' : -1.0278,
      'Jets01PassPtG10' : -1.0250,
      },
    },
  'w' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    },
  'z' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    },
  }
errors.JER = {
  'gg' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : -1.0166,
      'Jet2CutsGFPass' : -1.0126,
      'Jet2CutsVBFPass' : 1.0321,
      'Jets01FailPtG10' : -1.0004,
      'Jets01PassPtG10' : 1.0012,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : -1.0131,
      'Jet2CutsGFPass' : -1.0101,
      'Jet2CutsVBFPass' : -1.0273,
      'Jets01FailPtG10' : 1.0003,
      'Jets01PassPtG10' : 1.0009,
      },
    },
  'vbf' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : -1.0111,
      'Jet2CutsGFPass' : -1.0074,
      'Jet2CutsVBFPass' : -1.0074,
      'Jets01FailPtG10' : -1.0167,
      'Jets01PassPtG10' : 1.0049,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : -1.0081,
      'Jet2CutsGFPass' : -1.0063,
      'Jet2CutsVBFPass' : -1.0054,
      'Jets01FailPtG10' : -1.0162,
      'Jets01PassPtG10' : 1.0052,
      },
    },
  'w' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    },
  'z' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    },
  }
errors.PU = {
  'gg' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : 1.00,
      'Jet2CutsGFPass' : 1.00,
      'Jet2CutsVBFPass' : 1.00,
      'Jets01FailPtG10' : 1.00,
      'Jets01PassPtG10' : 1.00,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : 1.0112,
      'Jet2CutsGFPass' : 1.0103,
      'Jet2CutsVBFPass' : 1.0222,
      'Jets01FailPtG10' : 1.051,
      'Jets01PassPtG10' : 1.0032,
      },
    },
  'vbf' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : 1.00,
      'Jet2CutsGFPass' : 1.00,
      'Jet2CutsVBFPass' : 1.00,
      'Jets01FailPtG10' : 1.00,
      'Jets01PassPtG10' : 1.00,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : 1.0089,
      'Jet2CutsGFPass' : 1.0049,
      'Jet2CutsVBFPass' : 1.006,
      'Jets01FailPtG10' : 1.0207,
      'Jets01PassPtG10' : 1.0027,
      },
    }
  }
errors.PUID = {
  'gg' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : 1.0109,
      'Jet2CutsGFPass' : 1.0158,
      'Jet2CutsVBFPass' : 1.0368,
      'Jets01FailPtG10' : 1.0007,
      'Jets01PassPtG10' : 1.0034,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : 1.0142,
      'Jet2CutsGFPass' : 1.0174,
      'Jet2CutsVBFPass' : 1.0386,
      'Jets01FailPtG10' : 1.0022,
      'Jets01PassPtG10' : 1.0050,
      },
    },
  'vbf' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : 1.0134,
      'Jet2CutsGFPass' : 1.0160,
      'Jet2CutsVBFPass' : 1.0324,
      'Jets01FailPtG10' : 1.0115,
      'Jets01PassPtG10' : 1.0108,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : 1.0155,
      'Jet2CutsGFPass' : 1.0168,
      'Jet2CutsVBFPass' : 1.0328,
      'Jets01FailPtG10' : 1.0129,
      'Jets01PassPtG10' : 1.0128,
      },
    },
  'w' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    },
  'z' : {
    '7TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    '8TeV' : {
      'Jet2CutsFailVBFGF' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsVBFPass' : None,
      'Jets01FailPtG10' : None,
      'Jets01PassPtG10' : None,
      },
    },
  }

errors.QCDScale = {
  'gg' : {
    '8TeV' : {
      'Jets01PassPtG10' : -1.0119,
      'Jets01FailPtG10' : 1.0438,
      'Jet2CutsVBFPass' : -1.1585,
      'Jet2CutsGFPass' : 1.1550,
      'Jet2CutsFailVBFGF' : -1.0658,
      },
    '7TeV' : {
      'Jets01PassPtG10' : -1.0079,
      'Jets01FailPtG10' : -1.0281,
      'Jet2CutsVBFPass' : -1.2424,
      'Jet2CutsGFPass' : -1.1608,
      'Jet2CutsFailVBFGF' : 1.0615,
      },
    },
  'vbf' : {
    '8TeV' : {
      'Jets01PassPtG10' : -1.0139,
      'Jets01FailPtG10' : -1.1585,
      'Jet2CutsVBFPass' : 1.0313,
      'Jet2CutsGFPass' : 1.0302,
      'Jet2CutsFailVBFGF' : -1.0422,
      },
    '7TeV' : {
      'Jets01PassPtG10' : 1.0042,
      'Jets01FailPtG10' : 1.0552,
      'Jet2CutsVBFPass' : -1.0487,
      'Jet2CutsGFPass' : -1.0381,
      'Jet2CutsFailVBFGF' : 1.0654,
      },
    },
  'w' : {
    '8TeV' : {
      'Jets01PassPtG10' : None,
      'Jets01FailPtG10' : None,
      'Jet2CutsVBFPass' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsFailVBFGF' : None,
      },
    '7TeV' : {
      'Jets01PassPtG10' : None,
      'Jets01FailPtG10' : None,
      'Jet2CutsVBFPass' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsFailVBFGF' : None,
      },
    },
  'z' : {
    '8TeV' : {
      'Jets01PassPtG10' : None,
      'Jets01FailPtG10' : None,
      'Jet2CutsVBFPass' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsFailVBFGF' : None,
      },
    '7TeV' : {
      'Jets01PassPtG10' : None,
      'Jets01FailPtG10' : None,
      'Jet2CutsVBFPass' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsFailVBFGF' : None,
      },
    },
  }

errors.UE = {
  'gg' : {
    '8TeV' : {
      'Jets01PassPtG10' : -1.0619,
      'Jets01FailPtG10' : 1.1687,
      'Jet2CutsVBFPass' : 1.4228,
      'Jet2CutsGFPass' : 1.2317,
      'Jet2CutsFailVBFGF' : -1.0847,
      },
    '7TeV' : {
      'Jets01PassPtG10' : -1.0660,
      'Jets01FailPtG10' : 1.1588,
      'Jet2CutsVBFPass' : 1.5915,
      'Jet2CutsGFPass' : -1.1856,
      'Jet2CutsFailVBFGF' : -1.0698,
      },
    },
  'vbf' : {
    '8TeV' : {
      'Jets01PassPtG10' : -1.0199,
      'Jets01FailPtG10' : 1.0622,
      'Jet2CutsVBFPass' : -1.1030,
      'Jet2CutsGFPass' : -1.0396,
      'Jet2CutsFailVBFGF' : 1.0429,
      },
    '7TeV' : {
      'Jets01PassPtG10' : -1.0519,
      'Jets01FailPtG10' : 1.0642,
      'Jet2CutsVBFPass' : -1.1420,
      'Jet2CutsGFPass' : -1.1004,
      'Jet2CutsFailVBFGF' : 1.0762,
      },
    },
  'w' : {
    '8TeV' : {
      'Jets01PassPtG10' : None,
      'Jets01FailPtG10' : None,
      'Jet2CutsVBFPass' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsFailVBFGF' : None,
      },
    '7TeV' : {
      'Jets01PassPtG10' : None,
      'Jets01FailPtG10' : None,
      'Jet2CutsVBFPass' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsFailVBFGF' : None,
      },
    },
  'z' : {
    '8TeV' : {
      'Jets01PassPtG10' : None,
      'Jets01FailPtG10' : None,
      'Jet2CutsVBFPass' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsFailVBFGF' : None,
      },
    '7TeV' : {
      'Jets01PassPtG10' : None,
      'Jets01FailPtG10' : None,
      'Jet2CutsVBFPass' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsFailVBFGF' : None,
      },
    },
  }


errors.PDF = {
  'gg' : {
    '8TeV' : {
      'Jets01PassPtG10' : 1.0799,
      'Jets01FailPtG10' : 1.0844,
      'Jet2CutsVBFPass' : 1.0578,
      'Jet2CutsGFPass' : 1.0715,
      'Jet2CutsFailVBFGF' : 1.0809,
      },
    '7TeV' : {
      'Jets01PassPtG10' : 1.0830,
      'Jets01FailPtG10' : 1.0868,
      'Jet2CutsVBFPass' : 1.0593,
      'Jet2CutsGFPass' : 1.0743,
      'Jet2CutsFailVBFGF' : 1.0801,
      },
    },
  'vbf' : {
    '8TeV' : {
      'Jets01PassPtG10' : 1.0197,
      'Jets01FailPtG10' : 1.0186,
      'Jet2CutsVBFPass' : 1.0378,
      'Jet2CutsGFPass' : 1.0212,
      'Jet2CutsFailVBFGF' : 1.0196,
      },
    '7TeV' : {
      'Jets01PassPtG10' : 1.0220,
      'Jets01FailPtG10' : 1.0214,
      'Jet2CutsVBFPass' : 1.0448,
      'Jet2CutsGFPass' : 1.0299,
      'Jet2CutsFailVBFGF' : 1.0204,
      },
    },
  'w' : {
    '8TeV' : {
      'Jets01PassPtG10' : None,
      'Jets01FailPtG10' : None,
      'Jet2CutsVBFPass' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsFailVBFGF' : None,
      },
    '7TeV' : {
      'Jets01PassPtG10' : None,
      'Jets01FailPtG10' : None,
      'Jet2CutsVBFPass' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsFailVBFGF' : None,
      },
    },
  'z' : {
    '8TeV' : {
      'Jets01PassPtG10' : None,
      'Jets01FailPtG10' : None,
      'Jet2CutsVBFPass' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsFailVBFGF' : None,
      },
    '7TeV' : {
      'Jets01PassPtG10' : None,
      'Jets01FailPtG10' : None,
      'Jet2CutsVBFPass' : None,
      'Jet2CutsGFPass' : None,
      'Jet2CutsFailVBFGF' : None,
      },
    },
  }


######################################
######################################
######################################
######################################
######################################

# Now for a tables
for iTable in range(2):
  errorList = errors.contents
  nErrors = 4
  errorList = errorList[iTable*nErrors:min((iTable+1)*nErrors,len(errorList))]
  tableStr = ""
  extraCols = len(errorList)
  tableStr += r"\begin{tabular}{|c|l|c|"+"c|"*extraCols+r"} \hline" + "\n"
  tableStr += r"\multicolumn{3}{|c|}{} & \multicolumn{"+str(extraCols)+r"}{|c|}{Systematic Error (Relative Systematic Error)} \\ \hline" + "\n"
  #tableStr += r"\begin{tabular}{|c|l|c|"+"cc|"*extraCols+r"} \hline" + "\n"
  #tableStr += r"\multicolumn{3}{|c|}{} & \multicolumn{"+str(2*extraCols)+r"}{|c|}{Systematic Error (Relative Systematic Error)} \\ \hline" + "\n"
  tableStr += r"Sample & Category & Efficiency"
  for error in errorList:
    if error == "MCStat":
      error = "MC Statistics"
    elif error == "QCDScale":
      error = "QCD Scale"
    tableStr += r" & %s" % error
    #tableStr += r" & \multicolumn{2}{%s}" % error
  tableStr += r"\\ \hline"+ "\n"
  for energy in energies:
    for ds in datasets:
      dsLabel = "GF"
      if "vbf" in ds:
        dsLabel = "VBF"
      tableStr += "\multirow{"+str(len(categories))+"}{*}"
      tableStr += "{%s %s} \n" % (dsLabel,energy.replace("TeV"," TeV"))
      # now on to normal stuff
      for cat,label in zip(categories,labels):
        beginLineLen = len(tableStr)
        tableStr += " & "+label + " &"
        efficiency,efficiencyErr = effReader(energy,ds,cat,mass)
        tableStr += (r" %.2f\%%" % (efficiency*100.)) + " &"
        for error in errorList:
          tableStr += "\n    "
          if error == "MCStat":
            tableStr +=  r" $\pm$%.2f\%% (%.2f\%%) " % (efficiencyErr*100.,efficiencyErr/efficiency*100.)
            #tableStr +=  r" $\pm$%.2f\%% & (%.2f\%%) " % (efficiencyErr*100.,efficiencyErr/efficiency*100.)
          else:
            err = errors[error][ds][energy][cat]
            err = abs(err)-1.
            tableStr +=  r" $\pm$%.2f\%% (%.2f\%%) " % (err*efficiency*100.,err*100.)
            #tableStr +=  r" $\pm$%.2f\%% & (%.2f\%%) " % (err*efficiency*100.,err*100.)
          if error != errorList[-1]:
            tableStr += "&"
          tableStr += "    % "+error+"  "
        tableStr += "\n"+r"    \\ "+ "\n"
      tableStr = tableStr[:-1] + r" \hline "+ "\n"
  tableStr += r"\end{tabular}" + "\n"
  
  print
  print tableStr
  print
