Code bases For Recommener System
========================

Contributors: 

>    Jimmy Lin 
	

Overall Achievement
----------
1. One class support vector machine learns the job application model with 42.5% implementation Note that this model is seriously underfitted since the dimensionality extremely large, comparing to the training data.

Extract feature from description of Jobs
----------





Extract feature from description of Users
----------
Domain size of each variables: 
>   UserID 1
>   
>   WindowID 1
>   
>   State 115
>   
>   DegreeType 7
>   
>   Major 1287
>   
>   WorkHistoryCount 1
>   
>   TotalYearsExperience 1
>   
>   CurrentlyEmployed 3
>   
>   ManagedOthers 2
>   
>   ManagedHowMany 1

#How to improve?
1. Merge the majors: represent majors with major categories, maybe 7
   categories, like science, arts, engineering
2. Any better way to merge the major? to less categories may lose the
   information

##Vectorize text decription of job
Within 2000 texts with initially over 5000 unique words

1. Convert all terms into lowercase
2. Eliminate stop words
3. Remove punctuations
4. Only consider linguistic stemmes for each tokens
5. Remove all timestamps 
6. Remove all datestamps
8. Remove IP address and website link
9. Remove all phone numbers
10. Remove term consisting of pure digits
11. Check validity of a word
12. Remove tenses, plurals (get stems)

Current implementation curtail the size of keywords to be about 1700. 

####How to improve?
1. further remove non-content-bearing high-frequency and low-frequency words 
2. extract word phrase (consider more context, not just single word)
3. further shrink the keywords set

##Merge MAJOR attribute
1. recognize double major format
   - by notations of / and &
2. 

#TODO:
1. merge the major
2. apply large dataset on one-class learning
3. employ some other better one-class learning method
4. try some versabbi stuff: incorporate feature in the matrix factorization



