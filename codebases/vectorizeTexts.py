############################################################
##  FILENAME:   tf_idf.py    
##  VERSION:    1.0
##  SINCE:      2014-03-01
##  AUTHOR: 
##       Jimmy Lin (xl5224) - JimmyLin@utexas.edu  
##
############################################################
##  Edited by MacVim
##  Documentation auto-generated by Snippet 
############################################################

import nltk.corpus
from nltk.corpus import stopwords
from nltk.corpus import words
from nltk.text import TextCollection
from nltk.text import Text
from nltk.stem.lancaster import LancasterStemmer
from nltk.stem.isri import ISRIStemmer
import heapq, csv, sys
import string, re

WORDS = words.words()
stop = stopwords.words('english')

def createTextFromTxtFile (fn):
    '''    
    Create text object from external text file
    Input: 
        fn - name of the text file
    Output:
        text 
    '''
    f = open (fn, 'rU')
    raw = f.read()
    tokens = nltk.word_tokenize(raw)
    text = nltk.Text(tokens)
    return text

def getTextFromString (textString, removeStopwords=True):
    '''
    Get text from the input string
    Input:
        textString - text string
        removeStopwords - specify whether to remove stop words
    Output:
        text - Text object containing the input string
    '''
    tmpStr = textString.replace("\\r\\n", " ")
    tmpStr = tmpStr.replace("\\r", " ")
    tmpStr = tmpStr.replace("\\n", " ")
    tmpStr = tmpStr.replace("/", " ")
    tmpStr = tmpStr.replace("\\", " ") # 
    tmpStr = tmpStr.replace(".", " ") # period
    tmpStr = tmpStr.replace("'s", "") # possession
    tokens = nltk.word_tokenize(tmpStr)
    if removeStopwords:
        ptokens = [i for i in tokens if i not in stop]
    text = nltk.Text(ptokens)
    return text, ptokens

def getTextCollectionFromTxtFile (fn):
    '''
    Create text collection from external text files
    Input:
        fn - name of the external text file
    Output:
        textCollection containing all texts in the given file
    '''
    f = open (fn, 'rU')
    tc = []
    alltokens = []
    for line in f:
        text, tokens = getTextFromString(line)
        tc.append(text)
        alltokens.extend(tokens)
    return tc, alltokens

def preFilter(tokens):
    """
    Remove terms recognized as 
        timestamp
        datestamp
        phone number
        hyperlink
    """
    remainedTokens = []
    for t in tokens:
        if re.search(r'(.*):(.*):(.*)', t,  re.M|re.I): # time stamp
            continue
        elif re.search(r'(.*):(.*)(am|pm)', t,  re.M|re.I): # time stamp
            continue
        elif re.search(r'20(.*)-(.*)-(.*)', t,  re.M|re.I): # date stamp
            continue
        elif re.search(r'20(.*)/(.*)/(.*)', t,  re.M|re.I): # date stamp
            continue
        elif re.search(r'(.*)-(.*)-(.*)', t,  re.M|re.I): # phone number
            continue
        elif re.search(r'(.*)\.(.*)\.(.*)', t,  re.M|re.I): # hyperlink
            continue
        elif str.isdigit(t):
            continue
        elif str.isdigit(t.strip("-,")):
            continue
        elif len(t) <= 2: 
            continue
        elif t in stop:
            continue
        else:
            remainedTokens.append(t)
    return remainedTokens

def getStemmes (tokens):
    stemmedTokens = []
    st = LancasterStemmer()
    irsist = ISRIStemmer()
    for t in tokens:
        if t in WORDS: 
            stemmedTokens.append(st.stem(t))
    return stemmedTokens

def removePunctions (tokens):
    stripedTokens = []
    for t in tokens:
        tmp = t.strip(string.punctuation)
        tmp = tmp.replace("\xc2\xa0", " ") # non-break space
        tmp = tmp.replace("  ", " ") # redundant space
        tmps = tmp.split(' ')
        stripedTokens += tmps
    return stripedTokens

def removeOthers (tokens):
    stripedTokens = []
    for t in tokens:
        tmp = t
        tmp = tmp.strip(' ')
        tmps = tmp.split(' ')
        stripedTokens += tmps
    return stripedTokens

def lowerCase (tokens):
    caseIgnoredTokens = set([])
    for t in tokens:
        caseIgnoredTokens.add(t.lower())  
    return list(caseIgnoredTokens)

def unique (tokens):
    '''
    deduplicate the tokens
    '''
    nbefore = len(tokens)
    tokensAfterUniq = set(tokens)
    nafter = len(tokensAfterUniq)
    return list(tokensAfterUniq), nbefore, nafter

def computeTF (term, text, alltexts):
    '''
    Compute the TF from text collection
    Input:
        term - 
        text
        alltexts - text collections
    Output:
        tf
    '''
    return text.count(term)

class PriorityQueue:
    """
      Implements a priority queue data structure. Each inserted item
      has a priority associated with it and the client is usually interested
      in quick retrieval of the lowest-priority item in the queue. This
      data structure allows O(1) access to the lowest-priority item.

      Note that this PriorityQueue does not allow you to change the priority
      of an item.  However, you may insert the same item multiple times with
      different priorities.
    """
    def  __init__(self):
        self.heap = []
        self.count = 0

    def push(self, item, priority):
        entry = (priority, self.count, item)
        heapq.heappush(self.heap, entry)
        self.count += 1

    def pop(self):
        (_, _, item) = heapq.heappop(self.heap)
        return item

    def isEmpty(self):
        return len(self.heap) == 0


def prepruning (term): 
    """
    validate the word before counting its frequency in the documents
    """
    return False


def processTokens (tokens):
    lowerCasedTokens = lowerCase(tokens) 
    punctionRemovedTokens = removePunctions(lowerCasedTokens)
    stemmedTokens = getStemmes(punctionRemovedTokens)
    stampsRemovedTokens = preFilter(stemmedTokens)
    othersRemovedTokens = removeOthers(stampsRemovedTokens)
    return othersRemovedTokens

if __name__ == '__main__':
    """
    Configuration
    """
    inputTextsName = sys.argv[1]
    outputTextName = sys.argv[2]
    csvout = open(outputTextName, 'w+', 0)
    tfidfwriter = csv.writer(csvout, delimiter=' ')
    minSupport = 5 # min frequency for the word to be considered as a feature
    minDF = 10 # min number of documents the word needs to be in
    maxDFPercent = 0.90 # the max value of the expression (document frequency of a
    # word/total number of document) to be considered as good feature to be in
    # the document

    """
    Pre-processing
    """
    alltexts, alltokens = getTextCollectionFromTxtFile(inputTextsName)
    lowerCasedTokens = lowerCase(alltokens) 
    punctionRemovedTokens = removePunctions(lowerCasedTokens)
    stemmedTokens = getStemmes(punctionRemovedTokens)
    stampsRemovedTokens = preFilter(stemmedTokens)
    othersRemovedTokens = removeOthers(stampsRemovedTokens)
    uniqtokens, nbefore, nafter = unique(othersRemovedTokens)

    textCollection = TextCollection(alltexts)

    nTexts = len(alltexts)
    nUniqTokens = len(uniqtokens)
    print "nTexts: ", nTexts, "nUniqTokens: ", nUniqTokens

    allTF = []
    allTFIDF = []

    for ti in range(0, nTexts):
        text = alltexts[ti]
        text.tokens = processTokens (text.tokens) 
    
    for term in uniqtokens:
        tf_vector = []
        #tfidf_vector = []
        for ti in range(0, nTexts):
            text = alltexts[ti]
            tf = computeTF(term, text, textCollection)
            #tfidf = textCollection.tf_idf (term, text)
            tf_vector.append(tf)
            #tfidf_vector.append(tfidf)
        allTF.append(tf_vector)
        #allTFIDF.append(tfidf_vector)
        ## term processing
        support = sum(tf_vector)
        DF = len([i for i, e in enumerate(tf_vector) if e != 0])
        #DFpercent = max(tfidf_vector)
        #tfidfwriter.writerow([term, support, DF, DFpercent])
        tfidfwriter.writerow([term, support, DF])
        # post prunning
        if support < minSupport: continue
        else:
            print term, support, DF
        '''
        if DF < minDF: continue
        if DFpercent > maxDFPercent : continue
        '''
        ## pq.push((term, support, DF, DFpercent), -support)
    print "DONE"

