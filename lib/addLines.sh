#!/bin/bash
numLines=$(wc -l $1| cut -d " " -f 1)
let "numLines = -numLines/2"
nl -v $numLines $1 | awk '{print $2"\t"$1}' > $1.tmp