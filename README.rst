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

The main scripts you need to worry about is ``notlxbatch.sh``.  
It will write the 95% upper limit for each datacard to <category>.txt.out,
and the significance for each datacard to <category>.txt.sig

If you want to do combinations of categories, the ``combAllText.py`` script
makes the default combinations as .txt datacards.
*CombSplitAll* is the combination of all categories, *Jet2SplitCutsGFSplit* 
is the combination of the 2-jet categories,and *Jets01SplitCatAll* is the 
combination of the 0,1-jet categories.

