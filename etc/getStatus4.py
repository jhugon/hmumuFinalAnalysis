#!/usr/bin/python

import re
import os
import sys
import socket
import subprocess
import time
import datetime
import curses
import cStringIO

def findWindowCenter(window):
    yMaxWindow, xMaxWindow = window.getmaxyx()
    return yMaxWindow/2, xMaxWindow/2

def displayInWindowCenter(window,text,yShift=0):
    yCenter,xCenter = findWindowCenter(window)
    strLen = len(text)
    xShift = strLen/2
    for iCh, ch in enumerate(text):
      try:
          window.addch(yCenter-yShift,xCenter+iCh-xShift,ch)
      except:
          pass

def makeCenteredWindow(parent,ySize,xSize):
  yCenter,xCenter = findWindowCenter(parent)
  yBegin = yCenter-ySize/2
  xBegin = xCenter-xSize/2
  result = curses.newwin(ySize,xSize,yBegin,xBegin)
  return result

def drawBoxAroundWindow(parent,window,ch,spacing=2):
  yi,xi = window.getbegyx()
  yMaxWindow, xMaxWindow = window.getmaxyx()
  yf = yi + yMaxWindow
  xf = xi + xMaxWindow
  parent.vline(yi-spacing,xi-spacing,ch,yMaxWindow+2*spacing)
  parent.vline(yi-spacing,xf+spacing,ch,yMaxWindow+2*spacing)
  parent.hline(yi-spacing,xi-spacing,ch,xMaxWindow+2*spacing)
  parent.hline(yf+spacing,xi-spacing,ch,xMaxWindow+2*spacing+1)

def displayInWindowLeft(window,text,yShift=0,xShift=0,emph=False,emph2=False):
    yCenter,xCenter = findWindowCenter(window)
    for iCh, ch in enumerate(text):
      try:
          if emph2:
            window.addch(yCenter+yShift,iCh+xShift,ch,curses.color_pair(1) | curses.A_BOLD)
          elif emph:
            window.addch(yCenter+yShift,iCh+xShift,ch,curses.A_BOLD)
          else:
            window.addch(yCenter+yShift,iCh+xShift,ch)
      except:
          pass

def getJobsRunning():
  hostname = socket.gethostbyaddr(socket.gethostname())
  username = os.environ["LOGNAME"]
  njobs = 0
  njobsrun = 0
  if "cern" in hostname:
    p = subprocess.Popen(["bjobs"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output = p.communicate()[0]
    outputLines = output.splitlines()
    njobs = len(outputLines - 1)
    njobsrun = 0
    for line in outputLines:
      if "RUN" in line:
        njobsrun += 1
  else:  # Assume on HPC submit node
    p = subprocess.Popen(["qstat","-u",username,"-t"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output = p.communicate()[0]
    outputLines = output.splitlines()
    njobsQ  = 0
    for line in outputLines:
      if re.search(r"^[0-9].* Q ",line):
        njobsQ += 1
      if re.search(r"^[0-9].* R ",line):
        njobsrun += 1
    njobs = njobsrun + njobsQ
  return njobs,njobsrun

def curses_main(w):
    curses.init_pair(1,curses.COLOR_RED,curses.COLOR_BLACK)
    curses.init_pair(2,curses.COLOR_GREEN,curses.COLOR_BLACK)
    curses.init_pair(3,curses.COLOR_BLACK,curses.COLOR_RED)
    w2 = makeCenteredWindow(w,3,25)
    #w2.attrset(curses.color_pair(3))
    w2.attrset(curses.color_pair(2))
    w2.refresh()
    drawBoxAroundWindow(w,w2,"X")
    w.refresh()
    w3 = curses.newwin(1,27,0,0)
    noJobsSwitch = False
    try:
        while True:
            nSub, nRun = getJobsRunning()
            w2.erase()
            if nSub < 1:
              if not noJobsSwitch:
                curses.flash()
                noJobsSwitch = True
                w2.attrset(curses.color_pair(1) | curses.A_BOLD)
                w.attrset(curses.color_pair(3) | curses.A_BOLD)
                drawBoxAroundWindow(w,w2,"X")
                w.refresh()
            displayInWindowLeft(w2,"N Jobs:           "+str(nSub),-1)
            displayInWindowLeft(w2,"N Running:        "+str(nRun),0)
            if nSub > 0:
              displayInWindowLeft(w2,"Fraction Running: "+str(float(nRun)/nSub),1)
            else:
              displayInWindowLeft(w2,"No More Jobs!!",1)
            w2.refresh()
            displayInWindowLeft(w3,"As of: "+datetime.datetime.now().replace(microsecond=0).isoformat(' '),0)
            w3.refresh()
            time.sleep(300)
    finally:
        curses.nocbreak(); w.keypad(0); curses.echo()
        curses.curs_set(1)
        curses.endwin()

if __name__ == "__main__":
  if len(sys.argv) > 1:
    print getJobsRunning()
    sys.exit(0)
  try:
    curses.wrapper(curses_main)
  except KeyboardInterrupt:
    sys.exit(0)
