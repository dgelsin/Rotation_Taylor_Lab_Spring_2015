#!/usr/bin/env python
'''following this tutorial: http://synesthesiam.com/posts/an-introduction-to-pandas.html'''

import sys
import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab 
from pylab import show

f = sys.argv[1]
#import file from terminal (1 argument)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
data_table=DataFrame.from_csv(f, header=0, sep="\t")
#make a dataframe from imported file (this one is for tsv format, if I want  to get rid of the header put header=1)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
columns = data_table.columns
#make a list of the column names
#print columns
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
num_rows = len(data_table)
#count the number of rows
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#print data_table
#print the data frame
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#print columns
#print the column names
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#print num_rows
#print the # of rows in the dataframe
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#entire_column = data_table["Sample1"]
#print entire_column
#take an entire column from the dataframe and put it into a list, then print this list
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#two_entire_columns = data_table[["Sample1","Sample2"]]
#print two_entire_columns 
#do the same as just above, put put two columns into a list and print their values. You can do this with as many columns as long as you have the name of the column.
#You can get multiple columns out at the same time by passing in a list of strings
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#another_way_entire_column = data_table.Sample1
#print another_way_entire_column
#using dot syntax does that same thing as data_table["Sample1"]
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#five_items_in_column = data_table.Sample1.head()
#print five_items_in_column
# puts the first 5 rows in a column of a data frame into a list and prints it
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#n_number_of_rows_in_column = data_table.Sample1.head(56)
#print n_number_of_rows_in_column
#n_number_of_rows_in_column2 = data_table.Sample1.tail(66)
#print n_number_of_rows_in_column2
#putting a value in .head() gives you that specific number of rows. tail() is the same put gives you the last number of rows
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#lots_of_rows_in_column = data_table.Sample1.Sample2.Sample76.Sample115.head(79)
#print lots_of_rows_in_column
#can't do the dot syntax with lots of columns, but...
#lots_of_rows_in_column2 = data_table[["Sample1", "Sample2", "Sample76", "Sample115"]].head(79)
#print lots_of_rows_in_column2
#you can do it with this format :)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#data_table.columns = ["human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain", "human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain","human_brain", "chimp_brain", "rhesus_monkey_brain", "mouse_brain"]
#print data_table
#rename the column names, fuck yeah!
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#new_columns = data_table.human_brain.head(69)
#print new_columns
#and now you can take columns, put them into a list, using your own custom column names ;)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#old_columns = data_table.Sample76.head(4)
#print old_columns
#once you rename them, trying to pull out columns using old names won't work. You'd have the rename them back or get rid of the command that renamed them
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#std_dev_column = data_table.Sample76.std()
#print std_dev_column
#calculate standard deviation of a column
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
std_dev = data_table.std()
#print "Standard deviation for each column of the entire dataframe is " + str(std_dev)
#getting standard dev of each column
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#histogram = data_table.Sample1.hist()
#<matplotlib.axes.AxesSubplot at 0x4e54a50>
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#top = data_table.Sample1.values[88]
#print top
#inspect a single value from a column. The number in brackets specificies the row
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
empty = data_table.apply(lambda col: pd.isnull(col))
#print empty
#print empty.mC2.head(10)
#find if you have we have any missing values in the dataframe, specifically in a column with the later print statement
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
data_table.dropna()
#print data_table
#get rid of any rows that have NaN values
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
row7 = data_table.irow(7)
#print row7
#get a single row
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#num_of_no_expression = data_table[data_table.hM1 == 0]
#print num_of_no_expression
#returns every row in column 1 "hM1" with a value equal (==) 0 and the corresponding values in the same row in the other columns
#can also do <= or >= a value and select rows of interest that way
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
no_expression = data_table[data_table == 0]
#print no_expression
#this returns all columns in the datamatrix that have a value of 0 and every other value in the column = NaN
#no_expression = no_expression.fillna("")
#print no_expression
#returns same as one above, but replaces NaN with a blank space
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
non_expression = 0
#print data_table.apply(lambda r: r["hM1"] - r["hM2"], axis=1)
print data_table

#print data_table
#x = data_table
#plt.hist(x, range=None, bottom=None, histtype=u'bar', orientation='vertical',rwidth=None, log=False, color=None, label=None, stacked=False, hold=None)

