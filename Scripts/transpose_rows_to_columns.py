#!/usr/bin/env python
import sys

with open(sys.argv[1]) as f:
	lis = [x.split() for x in f]

for x in zip(*lis):
	for y in x:
		print (y + '\t')
		print ('\n')