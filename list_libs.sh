#!/bin/bash

for filename in `ls *.Rmd`;
do echo $filename;
   echo "";
   grep -h 'library(' $filename | sort | uniq;
   echo "";
done
