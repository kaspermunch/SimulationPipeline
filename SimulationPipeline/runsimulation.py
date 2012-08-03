#!/usr/bin/env python

import sys
sys.path.append( "/home/kmt/greatApes" )
sys.path.append( "/home/kmt/greatApes/scripts" )
sys.path.append( "/home/kmt/greatApes/CoalhmmPipeline_trunk")

from optparse import OptionParser
from SimulationPipeline import Utils

def main():

    usage="""%prog [options] [inputFile [outputFile]]

This program ....
It also does ..."""


    parser = OptionParser(usage=usage, version="%prog 1.0")

    parser.add_option("-v", "--verbose",
                      action="store_true",
                      dest="verbose",
                      #type="string",
                      default=False,
                      help="Print status output to STDERR")
    parser.add_option("--minprob",
#                      dest="minprob",
                      type="float",
                      default=0.5,
                      help="Print status output to STDERR")
    parser.add_option("--maxspan",
#                      dest="maxspan",
                      type="float",
                      default=float('inf'),
                      help="Print status output to STDERR")

    (options, args) = parser.parse_args()

    specFileName, pickleOutputFileName, coalhmmOptionsFile = args

    Utils.runSimulationsWithCoaSim(specFileName, pickleOutputFileName, options.minprob, options.maxspan, coalhmmOptionsFile)




if __name__ == "__main__":
    main()


