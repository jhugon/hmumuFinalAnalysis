#!/bin/bash

nice ./dataMC.py --n01Jet --energy 7TeV
nice ./dataMC.py --n01Jet --energy 8TeV
nice ./dataMC.py --n01JetMassOnly --energy 7TeV
nice ./dataMC.py --n01JetMassOnly --energy 8TeV
nice ./dataMC.py --n2Jet --energy 7TeV
nice ./dataMC.py --n2Jet --energy 8TeV
nice ./dataMC.py --n2JetMassOnly --energy 7TeV
nice ./dataMC.py --n2JetMassOnly --energy 8TeV
nice ./dataMC.py --n2JetVBFTight --energy 7TeV
nice ./dataMC.py --n2JetVBFTight --energy 8TeV
nice ./dataMC.py --n2JetVBFLoose --energy 7TeV
nice ./dataMC.py --n2JetVBFLoose --energy 8TeV
nice ./dataMC.py --n2JetGFTight --energy 7TeV
nice ./dataMC.py --n2JetGFTight --energy 8TeV

