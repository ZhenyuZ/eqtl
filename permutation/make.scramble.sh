#!/bin/bash
# run every bash command in scramble.cmd

while read -r line
do
	eval $line

done < "scramble.cmd"
