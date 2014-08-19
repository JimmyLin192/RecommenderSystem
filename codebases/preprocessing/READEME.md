Preprocessing Instruction
=======================

IMC experiment for winID=1
-------------------------------
1. generate matrix A

    python process_apps.py win1_Users.index job1.index ../../Dataset/apps.tsv A

2. generate matrix X (parameters are set within script)
    
    python process_users.py 

3. generate matrix Y
    
    python extractHTMLTxt.py ../../Dataset/splitjobs/jobs1.tsv
    python vectorizeTexts.py 
