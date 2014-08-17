#!/usr/bin/python

# run command: python $HOME/GBM_genomics/MAKE_VCF/split_fastq_threaded.py $pfx $splitSize
# $pfx is a barcode for the sample, and $splitSize is a command line argument

import os
import sys
import Queue
import threading
import time
from itertools import groupby, count
from Bio import SeqIO

def fastqSplit( fileName, pt, chunk=1000000 ):
    outFilePfx = "r" + str(pt)
    print fileName
    print outFilePfx
    with open( fileName, 'r' ) as datafile:
        groups = groupby( datafile, key=lambda k, line=count(): next(line) // chunk )
        n = 0
        for k, group in groups:
            with open( outFilePfx+"."+str(n), 'w' ) as outfile:
                outfile.write(''.join(group))
                outfile.close()
            n += 1
    datafile.close()

def removeUnpairedMates( sample, ptNum ):
    inFastQ = sample + "." + str(ptNum) + ".fastq"
    outFastQ = sample + "." + str(ptNum) + ".clean.fastq"
    print "Cleaning %s --> %s" % ( inFastQ, outFastQ )
    cleanFastq = open( outFastQ, 'w' )
    for rec in SeqIO.parse( open(inFastQ), 'fastq' ):
        mateInfo = rec.name[-2:]
        # print rec.name
        # print mateInfo
        if mateInfo == '/1' or mateInfo == '/2':
            SeqIO.write( rec, cleanFastq, 'fastq' )
    cleanFastq.close()
    # rename the "clean" file to the starting file
    os.rename( outFastQ, inFastQ )

# set up the global threading queue
threadQueue = Queue.Queue()
class FastQThread( threading.Thread ):

    def __init__( self, queue ):
        threading.Thread.__init__( self )
        self.queue = queue

    def run( self ):
        while True:
            # get the sample information from the queue
            info = self.queue.get()
            function = info[0]
            args = info[1:]
            try:
                function( *args )
            except:
                print "something wrong in " + self.name
                self.queue.task_done()
            # send the done signal
            self.queue.task_done()

startTime = time.time()
def main():

    # spawn two threads
    for i in range( 2 ):
        t = FastQThread( threadQueue )
        t.setDaemon( True )
        t.start()
    
    # first, clean up the files if the sizes don't match (indicative of unpaired reads)
    sampleBase = sys.argv[1]
    readsPerSplit = int(sys.argv[2])
    if os.path.getsize( sampleBase+".1.fastq" ) != os.path.getsize( sampleBase+".2.fastq" ):
        for i in [1, 2]:
            print (sampleBase+"."+str(i)+".fastq", sampleBase+"."+str(i)+".clean.fastq", readsPerSplit)
            threadQueue.put( (removeUnpairedMates, sampleBase, i) )
        # wait until everything is done
        print "waiting for cleaning to finish..."
        threadQueue.join()

    # send samples to the queue
    for i in [1, 2]:
        print (sampleBase+"."+str(i)+".fastq", i, readsPerSplit)
        threadQueue.put( (fastqSplit, sampleBase+"."+str(i)+".fastq", i, readsPerSplit) )
    # wait until everything is done
    print "waiting for splits to finish..."
    threadQueue.join()

## this is where we run the main function, main()
main()
print (time.time()-startTime)
sys.exit()
