# Author: Minwoo Bae
# Description: Preprocesses the job/user data sets
# Created 3/1/2014


import csv

from HTMLParser import HTMLParser
#from nltk.corpus import stopwords

#with open("stopwords/english") as stopword_file:
#	stopwords = stopword_file.readlines()
#	print stopwords

#stopwords = [line.strip() for line in open('stopwords/english')]
#print stopwords

with open("./../Dataset/sampled_jobs.tsv") as input:
	# print zip(*(line.strip().split('\t') for line in input))
	text = zip(*(line.strip().split('\t') for line in input))
	# print text
	# text = text[:, 3]

	#text2 = [row[3] for row in text]
	#print text2

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
#delete_list = ["\\r", "\\n", "\\r\\n", "\\t"]
output = ["" for n in range(len(text))] # is this OBOE???
for i in range (0, len(text)-1): # is this OBOE???
#	text[i].rstrip('\r')
#	text[i].rstrip('\n')

	output[i] = strip_tags(text[i])
	
	for word in delete_list:
		output[i] = output[i].replace(word, "")

	#output[i] = output[i].strip()
	#output[i] = output[i].lstrip('\n')
	#output[i] = output[i].rstrip('\n')
	#output[i] = output[i].strip('\n')
	#output[i] = output[i].rstrip('\r\n')

	# this doesn't work right now because each array element is a string, and w goes through characters - 
	#output[i] = [w for w in output[i] if not w in stopwords]

	print output[i]


output_file = open('output_file', 'w')
for item in output:
	output_file.write("%s\n" % item)


# defunct - 
# filtered_words = [w for w in word_list if not w in stopwords.words('english')]



