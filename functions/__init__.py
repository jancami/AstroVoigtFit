import os
from sys import platform

#this is the only environmental variable so far
if 'AVF_PYTHONDIR' in os.environ:
    if platform == "win32":
        PYTHONDIR = os.environ['AVF_PYTHONDIR']
    else:
        PYTHONDIR = os.environ['AVF_PYTHONDIR']
    
else:
    PYTHONDIR = os.path.dirname(__file__)
