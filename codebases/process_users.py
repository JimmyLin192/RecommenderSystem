#!/usr/bin/python2.7
############################################################
##    FILENAME:   process_users.py    
##    VERSION:    1.5
##    SINCE:      2014-03-24
##    AUTHOR: 
##        Jimmy Lin (xl5224) - JimmyLin@utexas.edu  
##
############################################################
##    Edited by MacVim
##    Documentation auto-generated by Snippet 
############################################################
import csv, sys
from util import *
from nltk.stem.lancaster import LancasterStemmer
from nltk.stem.isri import ISRIStemmer
from nltk.corpus import words
from vectorizeTexts import unique

WORDS = set(words.words())

def preprocessing (major):
    # preprocessing
    major = major.replace("&", " ")
    major = major.replace(".", " ")
    major = major.replace("/", " ")
    lowcase_major = major.lower()
    separated_majors = lowcase_major.split(" ")
    return separated_majors

def vectorize_major(major_dict):
    st = LancasterStemmer()
    major_list = []
    for major in major_dict.keys():
        separated_majors = preprocessing(major)
        for i in range(len(separated_majors)):
            if separated_majors[i] in WORDS:
                if len(separated_majors[i]) <= 3:
                    continue
                major_list.append(st.stem(separated_majors[i]))

    uniqtokens, nbefore, nafter = unique(major_list)

    final_dict = {}
    index = 0
    for major in uniqtokens:
        if final_dict.has_key(major):
            pass
        else:
            final_dict.update({major:index})
            index += 1
     
    return final_dict


with open('./../Dataset/users.tsv','rb') as tsvin, \
        open('./sampleData/newUsers.csv', 'wb') as csvTrainOut, \
        open('./process_users.log', 'wb') as log:

    original = sys.stdout
    sys.stdout = log
    tsvin = csv.reader(tsvin, delimiter='\t')
    csvTrainOut = csv.writer(csvTrainOut)

    """
    UserID,WindowID,Split,City,State,
    Country,ZipCode,DegreeType,Major,GraduationDate, 
    WorkHistoryCount,TotalYearsExperience,CurrentlyEmployed,ManagedOthers,ManagedHowMany
    """
    header = tsvin.next()
    nColumns = len(header)

    dictionaries = [] # dictionary for discretizing each column
    values = [] # numerical count of distinct label
    for i in range(0, nColumns):
        dictionaries.append({})
        values.append(0)

    ## DISCRITIZATION
    discretSet = range(3,9) + range(12,14) # set of index to discretize
    removeSet = [1,2,3,5,6,9] # set of index to remove
    remainedSet = list(set(range(0,len(header))).difference(set(removeSet)))
    remainedSet.sort()
    SCSUsers = [] # users matrix by simple coding scheme
    BCSUsers = [] # users matrix by binary coding scheme
    newheader = [header[i] for i in remainedSet]
    MAJOR_IDX = 8

    for user in tsvin:
        windowID = user[1]
        if not windowID == '1':
            break
        split = user[2]
        userid = user[0]
        for i in discretSet:
            ## restore the MAJOR value
            key = user[i]
            if not dictionaries[i].has_key(key):
                dictionaries[i][key] = values[i]
                values[i] += 1
            if not i == MAJOR_IDX:
                user[i] = dictionaries[i][key]

        scsUser = user
        #print scsUser
        SCSUsers.append(scsUser)

    """
    STEP: MERGE MAJOR
    """
    keywords_dict = vectorize_major(dictionaries[MAJOR_IDX])
    values[MAJOR_IDX] = len(keywords_dict.keys())
    nSCSFeatures = len(SCSUsers[0])
    print "nSCSFeatures:", nSCSFeatures
    for item in keywords_dict.items():
        print item
    
    print "Domain size of each variables: "
    newheader = []
    for i in range(0, nColumns):
        if i in removeSet: continue
        if values[i] <= 0: values[i] = 1
        print "  ", header[i], values[i]
        newheader = newheader + [header[i]] * values[i]

    """
    STEP: 
    """
    ## output the accumulated counter
    for i in range(nColumns):
        if i in removeSet: continue
        domainSize = len(dictionaries[i].items())
        if not (domainSize == 0):
            print '===========', header[i], '==========='
            items = dictionaries[i].items()
            if domainSize <= 300 or True:
                for key, value in items:
                    print key, ":::",  value

    ## HEADER FORMULATION
    sys.stdout = original
    task_size = len(SCSUsers)
    print task_size
    pb = ProgressBar(task_size)
    pb_index = 0
    
    csvTrainOut.writerows([newheader])
    nBCSFeatures = None
    st = LancasterStemmer()
    for scsUser in SCSUsers:
        bcsUser = [] # new user vector in binary coding scheme
        for i in range(nSCSFeatures):
            if i in removeSet: # ignore feature to be removed
                continue
            if i in discretSet: # 
                temp = [0 for x in range(0,values[i])]
                if i == MAJOR_IDX:
                    ## TODO: major representation
                    sep_majors = preprocessing(scsUser[MAJOR_IDX])
                    for sep_major in sep_majors:
                        word = st.stem(sep_major)
                        if keywords_dict.has_key(word):
                            temp[keywords_dict[word]] = 1
                else:
                    temp[scsUser[i]] = 1
                bcsUser = bcsUser + temp
            else: #
                bcsUser.append(scsUser[i])
        ## display the binary encoding result and guarantee its safety
        if nBCSFeatures is None:
            nBCSFeatures = len(bcsUser)
            print "nBCSFeatures (total):", nBCSFeatures
        else:
            assert(nBCSFeatures == len(bcsUser))
        ## output to external csv, piece by piece
        csvTrainOut.writerows([bcsUser])

        pb_index += 1
        pb.update(pb_index)
        pb.display()
