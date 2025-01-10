# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from interface import *

### This script allows to generate batch files that submit jobs with `htc-submit`

# Environment to load
activate = '. /afs/ifh.de/group/pitz/data/lixiangk/work/apps/python3/3.9.18/bin/activate'
#activate = 'module add python/3.9'

command = 'python'
script = 'TripletTuning.py'

job = CondorJob(command = command)

jobname = f'temp.sh'
args = [1] # if no input arguments, set it to "[]"
# if `submit = True`, the job will be submitted automatically, but this works only 
# by running the script from the `htc-submit` server
submit = False 

job.create(jobName = jobname, 
           inputName = script, 
           submit = submit, 
           args = args, 
           activate = activate)