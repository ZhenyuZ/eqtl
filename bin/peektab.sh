#!/bin/bash

head -n+3 $1 | awk -F"\t" -f ~/bin/transpose2.awk | cat -n
