#!/usr/bin/env python

import optparse
parser = optparse.OptionParser(description="Makes Limit Plots from output text from combine tool.")
parser.add_option("--xs", help="Limit on XS*BR",action="store_true",default=True)
args, fakeargs = parser.parse_args()

from helpers import *
import ROOT as root
import glob
import re
import os.path
from copy import deepcopy
from array import array
import sys

from xsec import *

##############################################
# Only use Matplotlib if not in CMSSW_6*
cmsswVersion = ""
mplGood = True
if os.environ.has_key("CMSSW_VERSION"):
  cmsswVersion = os.environ["CMSSW_VERSION"]
#def ComparePlotTable(*arg):
#  #Dummy so that there are no errors
#  pass
if "CMSSW_6" in cmsswVersion:
  mplGood = False
if mplGood:
  from etc.evilMatPlotlibFunctions import *
##############################################


def getData():
  #  [xNum,obs,low2sig,low1sig,median,high1sig,high2sig]

   # 1 sigma band
   oneSigBand = root.TGraphAsymmErrors()
   oneSigBand.SetName("")
   oneSigBand.SetTitle("")
   oneSigBand.SetFillColor(root.kGreen) 
   oneSigBand.SetLineStyle(2)
   oneSigBand.SetPoint(0,120,0.0437539)
   oneSigBand.SetPointError(0,0,0,0.01217779,0.01702032)
   oneSigBand.SetPoint(1,121,0.0425998)
   oneSigBand.SetPointError(1,0,0,0.01185657,0.01657137)
   oneSigBand.SetPoint(2,122,0.04123918)
   oneSigBand.SetPointError(2,0,0,0.01147788,0.01604209)
   oneSigBand.SetPoint(3,123,0.04017265)
   oneSigBand.SetPointError(3,0,0,0.01118104,0.01562721)
   oneSigBand.SetPoint(4,124,0.03897402)
   oneSigBand.SetPointError(4,0,0,0.01084743,0.01516094)
   oneSigBand.SetPoint(5,125,0.03786301)
   oneSigBand.SetPointError(5,0,0,0.01053821,0.01472876)
   oneSigBand.SetPoint(6,126,0.0365721)
   oneSigBand.SetPointError(6,0,0,0.01017892,0.01422659)
   oneSigBand.SetPoint(7,127,0.03568101)
   oneSigBand.SetPointError(7,0,0,0.009930907,0.01387996)
   oneSigBand.SetPoint(8,128,0.03426255)
   oneSigBand.SetPointError(8,0,0,0.009536114,0.01332817)
   oneSigBand.SetPoint(9,129,0.03325272)
   oneSigBand.SetPointError(9,0,0,0.009255053,0.01293535)
   oneSigBand.SetPoint(10,130,0.03198104)
   oneSigBand.SetPointError(10,0,0,0.008901113,0.01244066)
   oneSigBand.SetPoint(11,131,0.03126434)
   oneSigBand.SetPointError(11,0,0,0.008701637,0.01216186)
   oneSigBand.SetPoint(12,132,0.03021892)
   oneSigBand.SetPointError(12,0,0,0.00841067,0.01175519)
   oneSigBand.SetPoint(13,133,0.02926491)
   oneSigBand.SetPointError(13,0,0,0.008145146,0.01138408)
   oneSigBand.SetPoint(14,134,0.02841425)
   oneSigBand.SetPointError(14,0,0,0.007908387,0.01105318)
   oneSigBand.SetPoint(15,135,0.02770756)
   oneSigBand.SetPointError(15,0,0,0.007711698,0.01077827)
   oneSigBand.SetPoint(16,136,0.02687298)
   oneSigBand.SetPointError(16,0,0,0.007479413,0.01045362)
   oneSigBand.SetPoint(17,137,0.02620224)
   oneSigBand.SetPointError(17,0,0,0.00729273,0.0101927)
   oneSigBand.SetPoint(18,138,0.02564425)
   oneSigBand.SetPointError(18,0,0,0.007137427,0.009975642)
   oneSigBand.SetPoint(19,139,0.02503824)
   oneSigBand.SetPointError(19,0,0,0.00696876,0.009739904)
   oneSigBand.SetPoint(20,140,0.02459486)
   oneSigBand.SetPointError(20,0,0,0.006845355,0.009567427)
   oneSigBand.SetPoint(21,141,0.02405088)
   oneSigBand.SetPointError(21,0,0,0.006693953,0.00935582)
   oneSigBand.SetPoint(22,142,0.02357943)
   oneSigBand.SetPointError(22,0,0,0.006562736,0.009172424)
   oneSigBand.SetPoint(23,143,0.02327908)
   oneSigBand.SetPointError(23,0,0,0.006479141,0.009055588)
   oneSigBand.SetPoint(24,144,0.02285773)
   oneSigBand.SetPointError(24,0,0,0.00636187,0.008891683)
   oneSigBand.SetPoint(25,145,0.02262147)
   oneSigBand.SetPointError(25,0,0,0.006296114,0.008799779)
   oneSigBand.SetPoint(26,146,0.0223031)
   oneSigBand.SetPointError(26,0,0,0.006207502,0.008675931)
   oneSigBand.SetPoint(27,147,0.02211447)
   oneSigBand.SetPointError(27,0,0,0.006155003,0.008602555)
   oneSigBand.SetPoint(28,148,0.02201191)
   oneSigBand.SetPointError(28,0,0,0.006126457,0.008562658)
   oneSigBand.SetPoint(29,149,0.0217154)
   oneSigBand.SetPointError(29,0,0,0.006043931,0.008447316)
   oneSigBand.SetPoint(30,150,0.02164403)
   oneSigBand.SetPointError(30,0,0,0.006024068,0.008419554)

   # 2 sigma band
   twoSigBand = root.TGraphAsymmErrors()
   twoSigBand.SetName("")
   twoSigBand.SetTitle("")
   twoSigBand.SetFillColor(root.kYellow) 
   twoSigBand.SetLineStyle(2)
   twoSigBand.SetPoint(0,120,0.0437539)
   twoSigBand.SetPointError(0,0,0,0.02001524,0.03699274)
   twoSigBand.SetPoint(1,121,0.0425998)
   twoSigBand.SetPointError(1,0,0,0.0194873,0.03601697)
   twoSigBand.SetPoint(2,122,0.04123918)
   twoSigBand.SetPointError(2,0,0,0.01886489,0.03486661)
   twoSigBand.SetPoint(3,123,0.04017265)
   twoSigBand.SetPointError(3,0,0,0.018377,0.03396489)
   twoSigBand.SetPoint(4,124,0.03897402)
   twoSigBand.SetPointError(4,0,0,0.01782869,0.03295148)
   twoSigBand.SetPoint(5,125,0.03786301)
   twoSigBand.SetPointError(5,0,0,0.01732046,0.03201215)
   twoSigBand.SetPoint(6,126,0.0365721)
   twoSigBand.SetPointError(6,0,0,0.01672993,0.03092072)
   twoSigBand.SetPoint(7,127,0.03568101)
   twoSigBand.SetPointError(7,0,0,0.0163223,0.03016733)
   twoSigBand.SetPoint(8,128,0.03426255)
   twoSigBand.SetPointError(8,0,0,0.01567342,0.02896806)
   twoSigBand.SetPoint(9,129,0.03325272)
   twoSigBand.SetPointError(9,0,0,0.01521148,0.02811428)
   twoSigBand.SetPoint(10,130,0.03198104)
   twoSigBand.SetPointError(10,0,0,0.01462974,0.0270391)
   twoSigBand.SetPoint(11,131,0.03126434)
   twoSigBand.SetPointError(11,0,0,0.01430189,0.02643315)
   twoSigBand.SetPoint(12,132,0.03021892)
   twoSigBand.SetPointError(12,0,0,0.01382366,0.02554927)
   twoSigBand.SetPoint(13,133,0.02926491)
   twoSigBand.SetPointError(13,0,0,0.01338725,0.02474269)
   twoSigBand.SetPoint(14,134,0.02841425)
   twoSigBand.SetPointError(14,0,0,0.01299811,0.02402348)
   twoSigBand.SetPoint(15,135,0.02770756)
   twoSigBand.SetPointError(15,0,0,0.01267484,0.02342599)
   twoSigBand.SetPoint(16,136,0.02687298)
   twoSigBand.SetPointError(16,0,0,0.01229306,0.02272037)
   twoSigBand.SetPoint(17,137,0.02620224)
   twoSigBand.SetPointError(17,0,0,0.01198623,0.02215328)
   twoSigBand.SetPoint(18,138,0.02564425)
   twoSigBand.SetPointError(18,0,0,0.01173097,0.02168152)
   twoSigBand.SetPoint(19,139,0.02503824)
   twoSigBand.SetPointError(19,0,0,0.01145376,0.02116915)
   twoSigBand.SetPoint(20,140,0.02459486)
   twoSigBand.SetPointError(20,0,0,0.01125093,0.02079428)
   twoSigBand.SetPoint(21,141,0.02405088)
   twoSigBand.SetPointError(21,0,0,0.01100209,0.02033437)
   twoSigBand.SetPoint(22,142,0.02357943)
   twoSigBand.SetPointError(22,0,0,0.01078642,0.01993577)
   twoSigBand.SetPoint(23,143,0.02327908)
   twoSigBand.SetPointError(23,0,0,0.01064903,0.01968183)
   twoSigBand.SetPoint(24,144,0.02285773)
   twoSigBand.SetPointError(24,0,0,0.01045628,0.01932559)
   twoSigBand.SetPoint(25,145,0.02262147)
   twoSigBand.SetPointError(25,0,0,0.0103482,0.01912584)
   twoSigBand.SetPoint(26,146,0.0223031)
   twoSigBand.SetPointError(26,0,0,0.01020256,0.01885666)
   twoSigBand.SetPoint(27,147,0.02211447)
   twoSigBand.SetPointError(27,0,0,0.01011628,0.01869719)
   twoSigBand.SetPoint(28,148,0.02201191)
   twoSigBand.SetPointError(28,0,0,0.01006936,0.01861047)
   twoSigBand.SetPoint(29,149,0.0217154)
   twoSigBand.SetPointError(29,0,0,0.009933721,0.01835978)
   twoSigBand.SetPoint(30,150,0.02164403)
   twoSigBand.SetPointError(30,0,0,0.009901074,0.01829944)

   # Observed
   obsBand = root.TGraphAsymmErrors()
   obsBand.SetName("")
   obsBand.SetTitle("")
   obsBand.SetLineStyle(1)
   obsBand.SetPoint(0,120,0.04231358)
   obsBand.SetPointError(0,0,0,0,0)
   obsBand.SetPoint(1,121,0.0362299)
   obsBand.SetPointError(1,0,0,0,0)
   obsBand.SetPoint(2,122,0.03699423)
   obsBand.SetPointError(2,0,0,0,0)
   obsBand.SetPoint(3,123,0.0386155)
   obsBand.SetPointError(3,0,0,0,0)
   obsBand.SetPoint(4,124,0.03798886)
   obsBand.SetPointError(4,0,0,0,0)
   obsBand.SetPoint(5,125,0.03757517)
   obsBand.SetPointError(5,0,0,0,0)
   obsBand.SetPoint(6,126,0.03919741)
   obsBand.SetPointError(6,0,0,0,0)
   obsBand.SetPoint(7,127,0.03804189)
   obsBand.SetPointError(7,0,0,0,0)
   obsBand.SetPoint(8,128,0.03950008)
   obsBand.SetPointError(8,0,0,0,0)
   obsBand.SetPoint(9,129,0.02978723)
   obsBand.SetPointError(9,0,0,0,0)
   obsBand.SetPoint(10,130,0.02594321)
   obsBand.SetPointError(10,0,0,0,0)
   obsBand.SetPoint(11,131,0.02525321)
   obsBand.SetPointError(11,0,0,0,0)
   obsBand.SetPoint(12,132,0.02895005)
   obsBand.SetPointError(12,0,0,0,0)
   obsBand.SetPoint(13,133,0.04040421)
   obsBand.SetPointError(13,0,0,0,0)
   obsBand.SetPoint(14,134,0.04764681)
   obsBand.SetPointError(14,0,0,0,0)
   obsBand.SetPoint(15,135,0.04056238)
   obsBand.SetPointError(15,0,0,0,0)
   obsBand.SetPoint(16,136,0.02977736)
   obsBand.SetPointError(16,0,0,0,0)
   obsBand.SetPoint(17,137,0.02315553)
   obsBand.SetPointError(17,0,0,0,0)
   obsBand.SetPoint(18,138,0.01940371)
   obsBand.SetPointError(18,0,0,0,0)
   obsBand.SetPoint(19,139,0.02042459)
   obsBand.SetPointError(19,0,0,0,0)
   obsBand.SetPoint(20,140,0.02308895)
   obsBand.SetPointError(20,0,0,0,0)
   obsBand.SetPoint(21,141,0.02800417)
   obsBand.SetPointError(21,0,0,0,0)
   obsBand.SetPoint(22,142,0.02982162)
   obsBand.SetPointError(22,0,0,0,0)
   obsBand.SetPoint(23,143,0.03022439)
   obsBand.SetPointError(23,0,0,0,0)
   obsBand.SetPoint(24,144,0.03115675)
   obsBand.SetPointError(24,0,0,0,0)
   obsBand.SetPoint(25,145,0.03002883)
   obsBand.SetPointError(25,0,0,0,0)
   obsBand.SetPoint(26,146,0.02633442)
   obsBand.SetPointError(26,0,0,0,0)
   obsBand.SetPoint(27,147,0.02663513)
   obsBand.SetPointError(27,0,0,0,0)
   obsBand.SetPoint(28,148,0.02596744)
   obsBand.SetPointError(28,0,0,0,0)
   obsBand.SetPoint(29,149,0.02130271)
   obsBand.SetPointError(29,0,0,0,0)
   obsBand.SetPoint(30,150,0.02019224)
   obsBand.SetPointError(30,0,0,0,0)

   result = []
   X = root.Double()
   Y = root.Double()
   print "points: ",oneSigBand.GetN()
   for i in range(oneSigBand.GetN()):
     oneSigBand.GetPoint(i,X,Y)
     x = float(X)
     exp = float(Y)
     obsBand.GetPoint(i,X,Y)
     obs = float(Y)
     m1sig = exp-oneSigBand.GetErrorYlow(i)
     p1sig = exp+oneSigBand.GetErrorYhigh(i)
     m2sig = exp-twoSigBand.GetErrorYlow(i)
     p2sig = exp+twoSigBand.GetErrorYhigh(i)
     result.append([x,obs,m2sig,m1sig,exp,p1sig,p2sig])
   return result

class RelativePlot:
  def __init__(self,dataPoints, canvas, legend, caption, ylabel="95% CL Limit on #sigma/#sigma_{SM} (H#rightarrow#mu#mu)", xlabel="Integrated Luminosity [fb^{-1}]",caption2="",caption3="",ylimits=[],xlimits=[],vertLines=[],showObs=False,energyStr="8TeV",caption4=""):
    expGraph = root.TGraph()
    expGraph.SetLineStyle(2)
    expGraph.SetMarkerStyle(20)
    expGraph.SetMarkerSize(0.9)
    oneSigGraph = root.TGraphAsymmErrors()
    oneSigGraph.SetFillColor(root.kGreen)
    oneSigGraph.SetLineStyle(0)
    twoSigGraph = root.TGraphAsymmErrors()
    twoSigGraph.SetFillColor(root.kYellow)
    twoSigGraph.SetLineStyle(0)
    oneGraph = root.TGraph()
    oneGraph.SetLineColor(root.kRed)
    oneGraph.SetLineStyle(3)
    obsGraph = root.TGraph()
    obsGraph.SetLineColor(1)
    obsGraph.SetLineStyle(1)
    obsGraph.SetLineWidth(3)
    obsGraph.SetMarkerStyle(20)
    obsGraph.SetMarkerSize(1.1)
    obsGraph.SetMarkerColor(1)
    self.expGraph = expGraph
    self.oneSigGraph = oneSigGraph
    self.twoSigGraph = twoSigGraph
    self.oneGraph = oneGraph
    self.obsGraph = obsGraph
    iPoint = 0
    ymax = 1e-20
    ymin = 1e20
    for point in dataPoints:
      xNum = float(point[0])
      if len(xlimits)==2:
        if xNum<xlimits[0]:
            continue
        if xNum>xlimits[1]:
            continue
      obs = float(point[1])
      median = float(point[4])
      low2sig = median - float(point[2])
      low1sig = median - float(point[3])
      high1sig = float(point[5]) - median
      high2sig = float(point[6]) - median
      expGraph.SetPoint(iPoint,xNum,median) 
      oneSigGraph.SetPoint(iPoint,xNum,median) 
      twoSigGraph.SetPoint(iPoint,xNum,median) 
      oneSigGraph.SetPointError(iPoint,0.0,0.0,low1sig,high1sig) 
      twoSigGraph.SetPointError(iPoint,0.0,0.0,low2sig,high2sig) 
      oneGraph.SetPoint(iPoint,xNum,1.0)
      obsGraph.SetPoint(iPoint,xNum,obs)
      iPoint += 1

      if float(point[6]) > ymax:
        ymax = float(point[6])
      if obs > ymax:
        ymax = obs

    label = root.TLatex()
    #label.SetNDC()
    label.SetTextFont(root.gStyle.GetLabelFont("X"))
    label.SetTextSize(root.gStyle.GetLabelSize("X"))
    label.SetTextColor(root.kRed)
    label.SetTextAlign(22)
    self.label=label

    twoSigGraph.Draw("a3")
    twoSigGraph.GetXaxis().SetTitle(xlabel)
    twoSigGraph.GetYaxis().SetTitle(ylabel)
    twoSigGraph.GetYaxis().SetNdivisions(6,4,0)
    if len(ylimits)==2:
        twoSigGraph.GetYaxis().SetRangeUser(*ylimits)
    else:
        twoSigGraph.GetYaxis().SetRangeUser(0.0,ymax*1.1)
    if len(xlimits)==2:
        twoSigGraph.GetXaxis().SetRangeUser(*xlimits)
    oneSigGraph.Draw("3")
    expGraph.Draw("l")
    expGraph.Draw("p")
    if showObs:
      obsGraph.Draw("l")
      obsGraph.Draw("p")

    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    #tlatex.SetTextSize(0.05)
    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.88,caption2)
    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.82,caption3)
    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.88,caption4)

    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,caption)

    self.vertLine = root.TLine()
    self.vertLine.SetLineColor(root.kRed)
    self.vertLine.SetLineWidth(3)
    self.arrows = []
    for xPos in vertLines:
      self.vertLine.DrawLine(xPos,ylimits[0],xPos,ylimits[1])
      xAxis = twoSigGraph.GetXaxis()
      arrowLength = (xAxis.GetXmax() - xAxis.GetXmin())/10.
      arrowY = (ylimits[1]-ylimits[0])*0.8
      arrowHeadSize = 0.025
      arrow = root.TArrow(xPos,arrowY,xPos+arrowLength,arrowY,
                                arrowHeadSize,"|>")
      arrow.SetLineWidth(3)
      arrow.SetLineColor(root.kRed)
      arrow.SetFillColor(root.kRed)
      #arrow.SetAngle(40)
      arrow.Draw()
      self.arrows += [arrow]

    canvas.RedrawAxis()

if __name__ == "__main__":

  dirName = "statsInput/"
  
  outDir = "statsOutput/"
  
  root.gErrorIgnoreLevel = root.kWarning
  root.gROOT.SetBatch(True)
  setStyle()

  canvas = root.TCanvas()

  for period in ["8TeV"]:
    energyStr = "8TeV"
    energyStrWrite = energyStr
    if energyStr == "7P8TeV":
      energyStrWrite = "7 & 8 TeV"
    else:
      energyStrWrite = energyStr.replace("TeV"," TeV")
    caption2 = "#sqrt{s} = "+energyStrWrite
    legend = root.TLegend(0.58,0.70,0.9,0.9)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    for plotName in [""]:
      data = getData()
      xlimits = []
      ylimits = [0.,0.1]
      vertLines = []
      if len(data)<=1:
        continue
      xlabel="m_{H} [GeV/c^{2}]"
      ylabel="95% CL Limit on #sigma #times BR (H#rightarrow e^{+}e^{-}) [pb]"
      caption2 = "#sqrt{{s}} = 8 TeV L = {0:.1f} fb^{{-1}}".format(19.6)
      caption3 = ""
      caption4 = ""
      title = "Combination"
      incPlot = RelativePlot(data,canvas,legend,title,caption2=caption2,ylimits=ylimits,energyStr=energyStrWrite,xlabel=xlabel,caption3=caption3,showObs=True,xlimits=xlimits,vertLines = vertLines,ylabel=ylabel,caption4=caption4)
      saveAs(canvas,outDir+"xsbr_EE_"+energyStr)
