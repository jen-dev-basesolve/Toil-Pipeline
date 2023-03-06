import os
import tempfile

from toil.common import Toil
from toil.job import Job

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def gc_content(sequence):
    print()
    print('The GC content is..')

    
    seq=''
    for i in sequence:
        seq+=i.seq.rstrip()                                     # Remove trailing white space

    # Calculating GC content
    bac=round(gc_fraction(seq)*100,2)                           # gc_fraction is function imported to calculate GC content.
    return bac

def di_nuc(sequence):
    print()
    print('Dinucleotide frequencies..')
    
    
    seq=''
    for i in sequence:
        seq+=i.seq.rstrip()

    AA=round(100*seq.count_overlap('AA')/len(seq),2)            # count_overlap method counts AA in a sequence. overlap method count 3 for AA in 'AAA', whereas only .count method counts 1 AA for 'AAA'.
    TT=round(100*seq.count_overlap('TT')/len(seq),2)
    GG=round(100*seq.count_overlap('GG')/len(seq),2)
    CC=round(100*seq.count_overlap('CC')/len(seq),2)


    print('AA',AA,"%")
    print('TT',TT,"%")
    print('GG',GG,"%")
    print('CC',CC,"%")


if __name__ == "__main__":
    sequence = [i for i in SeqIO.parse('MonkeyPox.fasta','fasta')]
    jobstore: str = tempfile.mkdtemp("tutorial_quickstart") # tempfile.mkdtemp() makes a temporary directory.
                                                            # What does " jobstore: str = " do? 
    os.rmdir(jobstore)                                      # os.rmdir removes a directory. 
                                                            # rmdir() takes path as input, I am assuming jobstore from above line saves the path of the temporary file.
    options = Job.Runner.getDefaultOptions(jobstore)        # Job.Runner.getDefaultOptions takes in jobstore as input. 
                                                            # And returns options used by toil workflow as output in the form of argparse.Arguments 
    options.logLevel = "OFF"                                # when logLevel=OFF, only critical logs are shown
    options.clean = "always"                                # when clean = "always" the jobstore is deleted after completion. By default it is "never"

    hello_job = Job.wrapFn(gc_content,sequence)              # wrapFn is used to wrap the function, similar to decorators.
                                                            # Job.wrapFn(function,arguments).
                                                            # wrapFn has two arguments, 1st is the name of the function, 2nd is the inputs to that function
    second_job = Job.wrapFn(di_nuc,sequence)

    with Toil(options) as toil:                             # basically toil = Toil(options)
        print(toil.start(hello_job))                        # toil.start initiates toil workflow. Takes rootJob as input. Returns job value 

        print(toil.start(second_job))