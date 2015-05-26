#!/bin/bash

head -n+5 $1 | awk -F"\t" -f ~/bin/transpose2.awk | cat -n
