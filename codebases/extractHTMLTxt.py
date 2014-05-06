# Author: Jimmy Lin
# Description: Preprocesses the job/user data sets
# Created 3/1/2014


import csv
import sys
from HTMLParser import HTMLParser
#from nltk.corpus import stopwords

#with open("stopwords/english") as stopword_file:
#	stopwords = stopword_file.readlines()
#	print stopwords

#stopwords = [line.strip() for line in open('stopwords/english')]
#print stopwords

with open(sys.argv[1]) as input:
	text = zip(*(line.strip().split('\t') for line in input))

	text = text[3]
	text = text[1:len(text)-1]
	#print text

class MLStripper(HTMLParser):
	def __init__(self):
		self.reset()
		self.fed = []
	def handle_data(self, d):
		self.fed.append(d)
	def get_data(self):
		return ''.join(self.fed)

def strip_tags(html):
	s = MLStripper()
	s.feed(html)
	return s.get_data()

delete_list = ["\\r", "\\n", "\\r\\n", "\\t", "\\\r", "\\\n"]
output = ["" for n in range(len(text))] # is this OBOE???

for i in range (0, len(text)-1): # is this OBOE???
    try:
        output[i] = strip_tags(text[i])
    except:
        continue
	
	for word in delete_list:
		output[i] = output[i].replace(word, "")

	print output[i]

output_file = open('output_file', 'w')
for item in output:
	output_file.write("%s\n" % item)

