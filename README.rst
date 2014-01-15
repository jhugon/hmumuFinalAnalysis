H->mumu Final Analysis
======================

Analysis code for plotting distributions and final statistical analysis.
Requires ROOT ntuples created by hmumuAnalysis_.

.. _hmumuAnalysis: http://github.com/jhugon/hmumuAnalysis

This branch is an a move back to fitting dimuon mass from 110-160 GeV/c^2, 
possibly with a polynomial.

To use RooBernsteinFast implemination with order higher than 7, use this version of limit stuff:

::

  git clone https://github.com/jhugon/HiggsAnalysis-CombinedLimit.git
  git checkout higherOrderBernsteinFast

