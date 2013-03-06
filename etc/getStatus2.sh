#!/bin/bash

echo "==========================="

NJOBS=`ls *.txt | wc -l`
STARTTIME=`date +%s`

EXT=".out"
if [ -n "$1" ]; then
  EXT=$1
fi

echo "Looking for completed files *$EXT"

while true; do
  NCOMPLETE="0"
  FILELISTOUT=`ls *$EXT 2> /dev/null`
  if [ -n "$FILELISTOUT" ]; then
    for i in $FILELISTOUT; do
      NTMP=`cat $i | wc -l`
      if [ "$NTMP" -gt "0" ]; then
          NCOMPLETE=$(( $NCOMPLETE + 1 ))
      fi
    done
  fi
  NSTARTED=`ls -d Dir* 2> /dev/null | wc -l`
  echo "`date --rfc-3339=seconds` Jobs: $NJOBS Started: $NSTARTED Complete: $NCOMPLETE"
  if [ "$NCOMPLETE" -eq "$NJOBS" ]; then
    ENDTIME=`date +%s`
    echo "Took $(( $ENDTIME - $STARTTIME )) seconds"
    echo "All Jobs Complete"
    echo "==========================="
    exit "0"
  fi
  sleep 10
done

