#!/usr/bin/env python

import os,sys,inspect

def ensureImportInCurrentPath():
  currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
  parentdir = os.path.dirname(currentdir)
  sys.path.insert(0,currentdir) 
  sys.path.insert(0,parentdir) 

ensureImportInCurrentPath()
