Code bases For Recommener System
========================

Contributors: 

>    Jimmy Lin 
	
Overall Achievement
----------
1. [May 4th] One class support vector machine learns the job application model with **accuracy 42.5%**. Note that this model is seriously underfitted since the dimensionality extremely large, comparing to the training data.

2. [May 22nd] One class support vector machine learns the job applcation model
   with **accuracy 56.32%**. Note that in this implementation, we achieve
   **Major Vectorization** (1287 features) and **Job Description
   Vectorization** (2600 features). Training and testing with all application
   data in windowID=1 with about 300,000 training instances and 40,000 testing
   instances.

Extract feature from description of Jobs
----------
Domain size of each variables: 
> JobId 1
> 
> State 115
> 
> Description 1763
>

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

Extract feature from description of Users
----------
Domain size of each variables: 
> 
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
1. Any better way to merge the major? to less categories may lose the
   information

Current implementation curtail the size of keywords to be about 1700. 

##Merge MAJOR attribute
1. recognize double major format: by notations of / and &
2. vectorize each major by extracted keyword. The keyword extraction follows
   the same approach with that of user description
3. Other possible method:
    - **agglomerative clustering**, use vectorized data to measure proximity
      between majors and then cluster them. This method provides us very good
      scalability and gain any specified number of clusters. (attempt
      different cluster proximity measure.)
    
4. However, since no supervision, hard to evaluate the outcome of merging.
   One good solution is to **cross validate** on the various proximity measure.

#TODO:
1. employ some better one-class learning method
2. use uniform encoding scheme for state for both user and job
3. attempt versabbi learning






One-Class learning with different kernel
-------------------

1. In the truncated data, we found that radial kernel causes the best result.

2. 
