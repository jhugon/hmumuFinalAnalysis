#!/usr/bin/env python

import numpy
import matplotlib.pyplot as mpl

fig = mpl.figure()
ax = fig.add_subplot(111)

xRange = [0.5,15.0]
useCMSttbar = True

#From HATHOR: https://twiki.cern.ch/twiki/bin/view/Sandbox/CrossSectionsCalculationTool#Top_cross_sections_for_7_TeV
ttbarX = numpy.array([7,10,14])
ttbarY = numpy.array([148.5,375.6,833.1])
ttbarErr1 = numpy.array([18.2,45.0,97.0])
ttbarErr2 = numpy.array([10.0,18.0,27.0])
ttbarErr = numpy.sqrt(ttbarErr1**2 + ttbarErr2**2)

#From CMS
ttbarCMSX = numpy.array([7,8])
ttbarCMSY = numpy.array([157.5,225.197])

zX = numpy.array([7,8])
zY = numpy.array([3048,3503.71]) #m>50
#zY = numpy.array([2916,3351.97]) #60<m<120

if useCMSttbar:
  print("Using CMS ttbar Xsecs")
  ttbarX = ttbarCMSX
  ttbarY = ttbarCMSY
  mpl.plot(ttbarX,ttbarY,"g*")
else:
  mpl.errorbar(ttbarX,ttbarY,yerr=ttbarErr,linestyle="None",color="g")
mpl.plot(zX,zY,"bo")

ax.set_xlabel(r"$\sqrt{s}$ [TeV]")
ax.set_ylabel(r"$\sigma(pp \rightarrow t\bar{t}) [pb]$")
ax.set_title(r"$pp \rightarrow t\bar{t}$ Cross-Section v. $\sqrt{s}$")
ax.set_xlim(*xRange)

ttbarCoef,ttbarStuff = numpy.polynomial.polyfit(ttbarX,ttbarY,1,full=True)
zCoef,zStuff = numpy.polynomial.polyfit(zX,zY,1,full=True)
print ttbarStuff
print zStuff

ttbarModel = numpy.polynomial.Polynomial(ttbarCoef)
zModel = numpy.polynomial.Polynomial(zCoef)

xVals = numpy.linspace(*xRange)
ttbarPred = ttbarModel(xVals)
zPred = zModel(xVals)

mpl.plot(xVals,ttbarPred,"g-")
mpl.plot(xVals,zPred,"b-")

fig.savefig("PlotXS.png")

print("Z Coefs:     {}".format(zCoef))
print("ttbar Coefs: {}".format(ttbarCoef))

print("Z Xsec for 14 TeV:     {}".format(zModel(14.)))
print("ttbar Xsec for 14 TeV: {}".format(ttbarModel(14.)))

print("Z Xsec for 14 TeV (Scale Xsec):     {}".format(3503.71*14./8.))
print("ttbar Xsec for 14 TeV (Scale Xsec): {}".format(225.197*14./8.))
