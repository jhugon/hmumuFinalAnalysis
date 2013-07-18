#!/usr/bin/env python

# Inputs

# Local significance (in sigmas)
sigLocal = 3.26954

# Significance of reference point used to count 
# the number of of up crossings.
sig0 = 1.

# Number of up crossings counted at the above
# significance level
N0 = 2.


####################################################

from scipy import exp
from scipy.stats import chi2
from scipy.stats import norm

# Convert significance to p-value
pLocal = norm.sf(sigLocal)
p0 = norm.sf(sig0)

# Get the test statistic value corresponding to the p-value
u = chi2.isf(pLocal*2,1)
u0 = chi2.isf(p0*2,1)

# The main equations
N = N0 * exp(-(u-u0)/2.)
pGlobal = N + chi2.sf(u,1)/2.

# Further info
sigGlobal = norm.isf(pGlobal)
trialFactor = pGlobal/pLocal

print ("*"*70)
print("Input: sigLocal: %.2f sig0: %.2f N0: %.2f " % (sigLocal,sig0,N0))
print("local p-value corresponding to local significance: %.2e" % pLocal)
print("p-value corresponding to significance sig0: %.2e" % p0)
print
print("test stat values:")
print("  u: %.2f" % u)
print("  u0: %.2f" % u0)
print("  N: %.2f" % N)
print
print("global p-value:      %.5f" % pGlobal)
print("global significance: %.2f" % sigGlobal)
print("trial factor         %.2f" % trialFactor)
print ("*"*70)
