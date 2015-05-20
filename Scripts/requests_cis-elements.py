#!/usr/bin/env python

''''''''
#Requesting the YUE cis-linked gene webpage to find cis-regulatory element genomic coordinates that are linked to down- and up-regulated genes
''''''''
''''''''
#To find an API, right click a webpage in chrome and inspect element. 
#You will see the html format for the website. If you have an input box to 
#type in on the website (ie some where to type a gene name to get cis-elements linked to it)
#you need to open up each arrow (>) and click on the input box until it highlights something in the html code 
#This is your API code that you will feed into requests.get
''''''''

import sys
import requests

f=open(sys.argv[1])

genes = f.readlines()



gene_names = ['CACNG8', 'TUBB3', 'PACSIN1']

#print gene_names
#for i in genes:
#	fields = i.rstrip()
#	gene_names.append(fields)

for i in gene_names:
	r = requests.get('http://promoter.bx.psu.edu/ENCODE/get_human_linkage.php?assembly=hg19&gene='+i)
#the core of the request package. Basically the API url + whatever you want to feed into it (ie a gene name)

#print r.url
#print out how the url link will look like

#r.encoding = 'gb2312'
# print out the type of encoding your text file will be returned as. This can convert html format into rich text format, possible others?

print r.text
#this is how you return what is on a webpage

#print r.content


''''''
#reference API, this is what most API's will look like
''''''
#assembly=hg19&gene=
#http://promoter.bx.psu.edu/ENCODE/get_human_linkage.php?assembly=hg19&gene=SYNGR3

''''''
#feeding in parameters to url via dictionary
''''''
#with open(f, 'rU') as document:
	#dic = {} # initiate dictionary
	#for line in document:
		#key, value = line[:-1].split("\t") #split by tab delimiter
		#dic[key] = value

#should be in this format dic={'gene': 'SYNGR3'}

#you can also feed in a dictionary into the request.get function, pretty sweet!
#r = requests.get("http://promoter.bx.psu.edu/ENCODE/get_human_linkage.php?assembly=hg19&=", params=dic)


#http://docs.python-requests.org/en/latest/user/quickstart/ 
#request quickstart manual