import os
from sys import platform

#you will want to create some folders to keep all of your python projects in, i created a folder called home, 
#then within I created another folder called python, in that folder i store all of my python projects.


#set up environmental varibles, I used pycharm so its a little different, it varies system to system, I googled how to do it.
#create environmental variables,you need two: 
#PYTHONPATH = '/home/python'
#AVF_PYTHONDIR = '/home/python/AstroVoigtFit

#this is just a check.
if 'AVF_PYTHONDIR' in os.environ:
    if platform == "win32":
        PYTHONDIR = os.environ['AVF_PYTHONDIR']
    else:
        PYTHONDIR = os.environ['AVF_PYTHONDIR']
    
else:
    PYTHONDIR = os.path.dirname(__file__)
