import os
import tempfile
import requests
from toil.common import Toil
from toil.job import Job

# Goal is to pass a variable from indexable, which is a job for function fn().

def parent(job,kin):
    return kin


def raiseWrap(job,arg):
    var = arg + 112 
    job.log("".join(["=>>>>>>>>>>>>>>>>>>>> this is a variable from raiseWrap = ",str(var)]))
    return var

def fn(job,num):
    num = num + 12
    job.log("".join(["=>>>>>>>>>>>>>>>>>>>> this is num received from parent = ",str(num)]))
    return num

j1 = Job.wrapJobFn(parent,888)

returnedvalue = j1.rv()
j2 = j1.addChildJobFn(fn,returnedvalue)

j3 = j2.addChildJobFn(raiseWrap,j2.rv())


if __name__ == "__main__":
    jobstore: str = tempfile.mkdtemp("tutorial_promises")
    os.rmdir(jobstore)
    options = Job.Runner.getDefaultOptions(jobstore)
    options.logLevel = "INFO"

        
    with Toil(options) as toil:
        toil.start(j1)