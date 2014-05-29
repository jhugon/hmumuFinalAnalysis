H->mumu Final Analysis
======================

Analysis code for plotting distributions and final statistical analysis.
Requires ROOT ntuples created by hmumuAnalysis_.

.. _hmumuAnalysis: http://github.com/jhugon/hmumuAnalysis

This branch is an a move back to fitting dimuon mass from 110-160 GeV/c^2,  
with the former MSSM function: Exp(Breit-Wigner+1/m^2).  Background
function choice systematic is taken care of by adding an additive signal
systematic (not multiplicative like others).


Instructions
---------------

To get final limits, you must first create the statistics datacards, and 
then run the statistics program on them

To make statistics cards: ``./makeCards.py -m 125``, where the argument 
to ``-m`` is the higgs mass you want to assume.

``makeCards.py`` puts the statistics datacards in the ``statsCards/``
directory, along with scripts to help you do statistics studies.
Each category will have both a .txt and .root datacard file.  Both are 
required for the statistics software to work. 

The main scripts you need to worry about is ``findExpected.sh``.  
It will write the expected 95% upper limit for each datacard to <category>.txt.explimit,
and the expected significance for a 1*SM signal for each datacard to <category>.txt.expsig.
The expected significance output is straightforward, but the expected 95% upper limit output
shows the quantiles.  The 50% quantile is the median, so the number next to 50% is the *median
95% CLs upper limit*.

If you want to do combinations of categories, the ``combAllText.py`` script
makes the default combinations as .txt datacards.
*CombSplitAll* is the combination of all categories, *Jet2SplitCutsGFSplit* 
is the combination of the 2-jet categories,and *Jets01SplitCatAll* is the 
combination of the 0,1-jet categories.

For the limits, the results are all in terms of (sigma)/(sigma\_SM).  I usually ony report
down to the tenths place, e.g. the baseline limit for the 8 TeV 2-Jet combination for a Higgs
mass of 125 GeV is 10.2*SM.  With your new cuts, getting an 8 TeV 2-Jet combination
better than that is your goal.
