Code bases For Recommener System
========================

Contributors: 
    Jimmy Lin 
    Minwoo Bae
	

Overall Achievement
----------
1. One class support vector machine with 42.5% implementation
Note that this model is seriously underfitted since the dimensionality
extremely large, comparing to the training data.

Extract feature from decription of job
----------
Domain size of each variables: 
   UserID 1
   WindowID 1
   State 222
   DegreeType 7
   Major 47072
   WorkHistoryCount 1
   TotalYearsExperience 1
   CurrentlyEmployed 3
   ManagedOthers 2
   ManagedHowMany 1

#TODO:
1. Merge the majors: represent majors with major categories, maybe 7
   categories, like science, arts, engineering

#Vectorize text decription of job
Within 2000 texts with initially 3300 unique words

1. convert all terms into lowercase
2. eliminate stop words
3. Remove punctuations
4. only consider linguistic stemmes for each tokens
5. Remove all timestamps 
6. Remove all datestamps
8. remove IP and website link
9. Remove all phone numbers
10. Remove term consisting of pure digits
11. check validity of a word
12. remove tenses


#What we can improve?
1. non-content-bearing high-frequency and low-frequency words 
2. extract word phrase
