#!/bin/bash
find . -name "$1" |  xargs -I {} -i ksh -c 'grep -H -A 5 "Overall I/O System Response time distribution" {} | tail -1'
<<COMMENT
files=`find . -name "$1"` #| xargs -I {} -i ksh -c 'grep -H -A 8 "Overall I/O System Response time distribution" {} | tail -1'
#echo $files
for f in ${files[@]}
do
  BK=$IFS
  results=`grep -H -A 20 "Overall I/O System Response time distribution" $f`
  IFS=$'\n'
  for res in ${results[@]}
  do
    echo $res >temp
    echo $res
    #cut -d ' ' -f 2,5 temp
  done
  IFS=$BK
done
COMMENT

