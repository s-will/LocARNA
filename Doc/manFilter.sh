#!/bin/bash

file="$1"

ext="${file##*.}"

echo [TOC]
echo

echo "#" $(basename "$file" ".$ext")
echo

if [ "$ext" = "pl" ] ; then
  # if perl file than
  pod2markdown "$file"
else
  # directly convert man page
  # Note: grep and sed clean artifacts of the current conversion

  pandoc $file --from man --to markdown | grep -vx '\\' | sed 's/**\\-/**-/g'
fi
