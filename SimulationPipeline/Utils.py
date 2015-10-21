
import sys, os, glob
from SimulationPipeHotSpotCTMC import coaSimLeafHack
if coaSimLeafHack:
    sys.path.insert(0, '/Users/kasper/Desktop/coasim_trunk/Python/build/lib.macosx-10.5-x86_64-2.7')
import CoaSim
from CoaSim.popStructure import Population as P, Sample as S, Merge as M, Migration as Mi, Growth
from SimulationPipeHotSpotCTMC import *
from newick.tree import Leaf
from CoalhmmPipeline import Table
from MultiPurpose import Needle
from optparse import OptionParser

def runSimulationsWithCoaSimAndILS09Script():

    usage="""%prog [options] specFileName, pickleOutputFileName, coalhmmOptionsFile

This program runs the simulation pipeline and takes the following arguments:
 - a specification python file
 - the name of the file to write the pickled result to
 - an option file for coalhmm
 - an option file for bppseqgen
"""

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
                      help="Minimal probability flanking recombination window")
    parser.add_option("--maxspan",
#                      dest="maxspan",
                      type="float",
                      default=float('inf'),
                      help="Maximal size of recombination windows reported")

    (options, args) = parser.parse_args()

    if len(args) > 4:
        parser.error("Too many arguments")
    if len(args) < 4:
        parser.error("Too few arguments")

    specFileName, pickleOutputFileName, coalhmmOptionsFile, bppseqgenOptionsFile = args

    runSimulationsWithCoaSimAndILS09(specFileName, pickleOutputFileName, options.minprob, options.maxspan, coalhmmOptionsFile, bppseqgenOptionsFile)


def runSimulationsWithCoaSimAndILS16Script():

    usage="""%prog [options] specFileName, pickleOutputFileName, coalhmmOptionsFile

This program runs the simulation pipeline and takes the following arguments:
 - a specification python file
 - the name of the file to write the pickled result to
 - an option file for coalhmm
 - an option file for bppseqgen
"""

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
                      help="Minimal probability flanking recombination window")
    parser.add_option("--maxspan",
#                      dest="maxspan",
                      type="float",
                      default=float('inf'),
                      help="Maximal size of recombination windows reported")

    (options, args) = parser.parse_args()

    if len(args) > 4:
        parser.error("Too many arguments")
    if len(args) < 4:
        parser.error("Too few arguments")

    specFileName, pickleOutputFileName, coalhmmOptionsFile, bppseqgenOptionsFile = args

    runSimulationsWithCoaSimAndILS09(specFileName, pickleOutputFileName, options.minprob, options.maxspan, coalhmmOptionsFile, bppseqgenOptionsFile)



def runSimulationsWithCoaSimAndILSCTMCScript():

    usage="""%prog [options] specFileName, pickleOutputFileName, coalhmmOptionsFile

This program runs the simulation pipeline and takes the following arguments:
 - a specification python file
 - the name of the file to write the pickled result to
 - an option file for coalhmm
 - an option file for bppseqgen
"""

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
                      help="Minimal probability flanking recombination window")
    parser.add_option("--maxspan",
#                      dest="maxspan",
                      type="float",
                      default=float('inf'),
                      help="Maximal size of recombination windows reported")

    (options, args) = parser.parse_args()

    if len(args) > 3:
        parser.error("Too many arguments")
    if len(args) < 3:
        parser.error("Too few arguments")

    specFileName, pickleOutputFileName, bppseqgenOptionsFile = args

    runSimulationsWithCoaSimAndILSCTMC(specFileName, pickleOutputFileName, bppseqgenOptionsFile)




def runSimulationsWithMaCS(inp, outp, minProb, maxSpan, coalhmmOptionsFile, bppseqgenOptionsFile):
    """
    Simulate a number of 
    """

    assert os.path.exists(os.path.abspath(coalhmmOptionsFile)) and os.path.exists(os.path.abspath(bppseqgenOptionsFile)), "%s %s" % (coalhmmOptionsFile, bppseqgenOptionsFile)

    with open(os.path.abspath(coalhmmOptionsFile)) as f:
        coalhmmModelLines = sorted([" ".join(x.split()) for x in f.readlines() if x.startswith('model ') or x.startswith('rate_distribution ')])
    with open(os.path.abspath(bppseqgenOptionsFile)) as f:
        bppseqgenModelLines = sorted([" ".join(x.split()) for x in f.readlines() if x.startswith('model ') or x.startswith('rate_distribution ')])
    assert coalhmmModelLines == bppseqgenModelLines, "model and rate distributions lines not identical in coalhmm and bppseqgen option files"

    import imp
    SimSpec = imp.new_module("SimSpec")
    exec open(inp).read() in SimSpec.__dict__
    sys.modules[SimSpec] = "SimSpec"

    # simulate:
    simHook = MaCSSimulationHook()
    if "recMap" in dir(SimSpec):
        seq = simulateMaCS(length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u, addOutGroup=SimSpec.addOutGroup, \
                           T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123, \
                           N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123, \
                           optionsFile=os.path.abspath(bppseqgenOptionsFile), hook=simHook, recMap=SimSpec.recMap)
    else:
        seq = simulateMaCS(length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u, addOutGroup=SimSpec.addOutGroup, \
                           T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123, \
                           N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123, \
                           optionsFile=os.path.abspath(bppseqgenOptionsFile), hook=simHook)

    assert len(simHook.recombinationPoints) == len(simHook.trees) - 1, "%s vs %s" % (len(simHook.recombinationPoints), len(simHook.trees) - 1)

    # mapping between externa leaf label and state corresponding to the tree. note that
    # the simulated states do not distinguish state 0 and 1. so some 0s are really 1s
    stateMap = {"0": 0, "1": 3, "2": 2}

    # filter the recombination events to extract recombinations resulting in diffent topologies:
    recPoints = list()
    fromState = list()
    toState = list()
    simStates = list()

    def getTreeExternalLeaf(tree):

#         left, right = tree.get_edges()
#         if isinstance(left[0], Leaf):
#             return left[0].identifier
#         else:
#             return right[0].identifier

        left, right = tree.get_edges()
        if isinstance(left[0], Leaf) and left[0].identifier == '3': # mayby identifier is not an int but a str...
            # outgroup is included - call on ingroup
            return getTreeExternalLeaf(right[0])
        elif isinstance(right[0], Leaf) and right[0].identifier == '3': # mayby identifier is not an int but a str...
            # outgroup is included - call on ingroup
            return getTreeExternalLeaf(left[0])        
        elif isinstance(left[0], Leaf):
            return left[0].identifier
        else:
            return right[0].identifier

    for i in range(len(simHook.trees)-1):

        #######################################################
        leftCoalTimes = getCoalTimes(simHook.trees[i])
        rightCoalTimes = getCoalTimes(simHook.trees[i+1])
        t12 = SimSpec.T12 / SimSpec.g / (2*SimSpec.NeRef)
        leftState = stateMap[leftExtLeaf]
        if leftState == 0 and min(leftCoalTimes) > t12:
            leftState = 1
        rightState = stateMap[rightExtLeaf]
        if rightState == 0 and min(rightCoalTimes) > t12:
            rightState = 1
        if leftState != rightState:
            recPoints.append(simHook.recombinationPoints[i])
            recTimes.append(simHook.recombinationTimes[i])
            fromState.append(leftState)
            toState.append(rightState)
        #######################################################

#         leftExtLeaf, rightExtLeaf = getTreeExternalLeaf(simHook.trees[i]),  getTreeExternalLeaf(simHook.trees[i+1])
# #         print 'left', simHook.trees[i], leftExtLeaf
# #         print 'right', simHook.trees[i+1], rightExtLeaf
#         if leftExtLeaf != rightExtLeaf:
#             recPoints.append(simHook.recombinationPoints[i])
#             fromState.append(stateMap[leftExtLeaf])
#             toState.append(stateMap[rightExtLeaf])

    # run coalhmm on simulated sequences:	
    estimHook = ILS09estimationHook(minProb, maxSpan)
    t = estimate_ils09(seq, length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u,
                       addOutGroup=SimSpec.addOutGroup, T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123,
                       N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123,
                       optionsFile=os.path.abspath(coalhmmOptionsFile), hook=estimHook)

    stats = { 'simRecPoints': recPoints,
              'simFromState': fromState,
              'simToState': toState,
              'infRecPoints': [estimHook.startCoord[i] + (estimHook.endCoord[i] - estimHook.startCoord[i])/2.0 for i in range(len(estimHook.startCoord))],
              'infStartCoords': estimHook.startCoord,
              'infEndCoords': estimHook.endCoord,
              'infFromState': estimHook.fromState,
              'infToState': estimHook.toState,
              'simParameters': dict([(k, str(v[0])) for k, v in t.data.items() if k.startswith('_')]),
              'infParameters': dict([(k, v[0].strip()) for k, v in t.data.items() if not k.startswith('_')]),
              'seqLength': SimSpec.L }

    with open(outp, 'w') as f:
        pickle.dump(stats, f)

#    os.system("cp " + seq + "* tmp_macs/.")

    for f in glob.glob(seq + "*"):
        os.unlink(f)


def runSimulationsWithCoaSimAndILS09(inp, outp, minProb, maxSpan, coalhmmOptionsFile, bppseqgenOptionsFile):
    """
    Simulate a number of 
    """

    assert os.path.exists(os.path.abspath(coalhmmOptionsFile)) and os.path.exists(os.path.abspath(bppseqgenOptionsFile)), "%s %s" % (coalhmmOptionsFile, bppseqgenOptionsFile)

    with open(os.path.abspath(coalhmmOptionsFile)) as f:
        coalhmmModelLines = sorted([" ".join(x.split()) for x in f.readlines() if x.startswith('model ') or x.startswith('rate_distribution ')])
    with open(os.path.abspath(bppseqgenOptionsFile)) as f:
        bppseqgenModelLines = sorted([" ".join(x.split()) for x in f.readlines() if x.startswith('model ') or x.startswith('rate_distribution ')])
    assert coalhmmModelLines == bppseqgenModelLines, "model and rate distributions lines not identical in coalhmm and bppseqgen option files"

    import imp
    SimSpec = imp.new_module("SimSpec")
    exec open(inp).read() in SimSpec.__dict__
    sys.modules[SimSpec] = "SimSpec"
    
    def spec(**args):
        c1 = args["N1"] / args["NeRef"]
        c2 = args["N2"] / args["NeRef"]
        c3 = args["N3"] / args["NeRef"]
        c12 = args["N12"] / args["NeRef"]
        c123 = args["N123"] / args["NeRef"]
        t1 = args["T1"] / args["g"] / (2*args["NeRef"])
        t12 = args["T12"] / args["g"] / (2*args["NeRef"])
        t123 = args["T123"] / args["g"] / (2*args["NeRef"])
        speciesTree = P(c123, M(t12, [P(c3, S(1)), P(c12, M(t1, [P(c1, S(1)), P(c2, S(1))]))]))
        return speciesTree

    # simulate:
    simHook = CoaSimSimulationHook()
    if "recMap" in dir(SimSpec):
        seq = simulateCoasim(spec, length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u, addOutGroup=SimSpec.addOutGroup, \
                             T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123, \
                             N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123, \
                             optionsFile=os.path.abspath(bppseqgenOptionsFile), hook=simHook, recMap=SimSpec.recMap)
    else:
        seq = simulateCoasim(spec, length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u, addOutGroup=SimSpec.addOutGroup, \
                             T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123, \
                             N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123, \
                             optionsFile=os.path.abspath(bppseqgenOptionsFile), hook=simHook)

    assert len(simHook.recombinationPoints) == len(simHook.recombinationTimes), "%s vs %s" % (len(simHook.recombinationPoints), len(simHook.recombinationTimes))
    assert len(simHook.recombinationPoints) == len(simHook.trees) - 1, "%s vs %s" % (len(simHook.recombinationPoints), len(simHook.trees) - 1)


    # mapping between externa leaf label and state corresponding to the tree. note that
    # the simulated states do not distinguish state 0 and 1. so some 0s are really 1s
    stateMap = {"0": 0, "1": 3, "2": 2}

    # filter the recombination events to extract recombinations resulting in diffent topologies:
    recPoints = list()
    recTimes = list()
    recLeaves = list()
    fromState = list()
    toState = list()
    simStates = list()
    isSameTree = list()
    isSameTopology = list()
    isSameState = list()
    allRecPoints = list()
    allRecTimes = list()
    allRecLeaves = list()
    #simCoalTimes = list()

    def getTreeExternalLeaf(tree):
        left, right = tree.get_edges()
        if isinstance(left[0], Leaf):
            return left[0].identifier
        else:
            return right[0].identifier

    def getCoalTimes(tree):        
        th = TreeHeight()
        tree.dfs_traverse(th)
        return sorted(list(set(th.max - x for x in th.depths)))


    for i in range(len(simHook.trees)-1):
        leftExtLeaf, rightExtLeaf = getTreeExternalLeaf(simHook.trees[i]),  getTreeExternalLeaf(simHook.trees[i+1])

        #######################################################
        leftCoalTimes = getCoalTimes(simHook.trees[i])
        rightCoalTimes = getCoalTimes(simHook.trees[i+1])
        #simCoalTimes.append([x * SimSpec.g * 2 * SimSpec.NeRef for x in sorted(leftCoalTimes)])
        t12 = SimSpec.T12 / SimSpec.g / (2*SimSpec.NeRef)
        leftState = stateMap[leftExtLeaf]
        if leftState == 0 and min(leftCoalTimes) > t12:
            leftState = 1
        rightState = stateMap[rightExtLeaf]
        if rightState == 0 and min(rightCoalTimes) > t12:
            rightState = 1
        if leftState != rightState:
            recPoints.append(simHook.recombinationPoints[i])
            recTimes.append(simHook.recombinationTimes[i])
            if coaSimLeafHack:
                # commented out because this functionality has a memory leak in CoaSim
                recLeaves.append(simHook.recombinationLeaves[i])
            fromState.append(leftState)
            toState.append(rightState)
        #######################################################


        isSameTree.append(int(str(simHook.trees[i]) == str(simHook.trees[i+1])))
        isSameTopology.append(int(leftExtLeaf == rightExtLeaf))
        isSameState.append(int(leftState == rightState))
        allRecPoints.append(simHook.recombinationPoints[i])
        allRecTimes.append(simHook.recombinationTimes[i])
        if coaSimLeafHack:
            # commented out because this functionality has a memory leak in CoaSim
            allRecLeaves.append(simHook.recombinationLeaves[i])


    # run coalhmm on simulated sequences:	
    estimHook = ILS09estimationHook(minProb, maxSpan)
    t = estimate_ils09(seq, length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u,
                       addOutGroup=SimSpec.addOutGroup, T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123,
                       N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123,
                       optionsFile=os.path.abspath(coalhmmOptionsFile), hook=estimHook)

    ## estimHook = CTMCestimationHook()
    ## t = estimate_3sp_ctmc(seq, r=r, g=g, u=u, T1=T1, T2=T12, N1=N1, nstates=10, hook=estimHook)
    ## t = estimate_untitled_3_all(seq, T1 * u, T12 * u, 1.0/(N1 * u * 2*g), r/(u*g), nstates=10)
    ## t = estimate_untitled_3_all(seq, T1 * u, T12 * u, 1.0/(N1 * u * 2*g), r/(u*g), nstates=20, suffix= "_%s_%d" % (tag, iteration))
    ## os.system("cp %s sequence1.fasta" % seq)

    if not coaSimLeafHack:
        recLeaves = [-1] * len(recTimes)
        allRecLeaves = [-1] * len(recLeaves)

    if t is not None:
        stats = { 'ilsBases': estimHook.ilsBases,
                  'non_ilsBases': estimHook.non_ilsBases,
                  'simRecPoints': recPoints,
                  'simRecTimes': recTimes,
                  'simRecLeaves': recLeaves,
                  'simFromState': fromState,
                  'simToState': toState,
                  #'simCoalTimes': simCoalTimes,
                  'infRecPoints': [estimHook.startCoord[i] + (estimHook.endCoord[i] - estimHook.startCoord[i])/2.0 for i in range(len(estimHook.startCoord))],
                  'infStartCoords': estimHook.startCoord,
                  'infEndCoords': estimHook.endCoord,
                  'infFromState': estimHook.fromState,
                  'infToState': estimHook.toState,
                  'allSimRecPoints': allRecPoints,
                  'allSimRecTimes': allRecTimes,
                  'allSimRecLeaves': allRecLeaves,
                  'isSameTree': isSameTree,
                  'isSameTopology': isSameTopology,
                  'isSameState': isSameState,
                  'simParameters': dict([(k, str(v[0])) for k, v in t.data.items() if k.startswith('_')]),
                  'infParameters': dict([(k, v[0].strip()) for k, v in t.data.items() if not k.startswith('_')]),
                  'seqLength': SimSpec.L }
    #               'simParameters': dict([(k, v) for k, v in t.data.items() if k.startswith('_')]),
    #               'infParameters': dict([(k, v) for k, v in t.data.items() if not k.startswith('_')])}

        with open(outp, 'w') as f:
            pickle.dump(stats, f)

#         os.system("cp " + seq + "* .")
#         print seq
#         sys.exit()

    for f in glob.glob(seq + "*"):
        os.unlink(f)


def runSimulationsWithCoaSimAndILS16(inp, outp, minProb, maxSpan, coalhmmOptionsFile, bppseqgenOptionsFile):
    """
    Simulate a number of
    """

    assert os.path.exists(os.path.abspath(coalhmmOptionsFile)) and os.path.exists(os.path.abspath(bppseqgenOptionsFile)), "%s %s" % (coalhmmOptionsFile, bppseqgenOptionsFile)

    with open(os.path.abspath(coalhmmOptionsFile)) as f:
        coalhmmModelLines = sorted([" ".join(x.split()) for x in f.readlines() if x.startswith('model ') or x.startswith('rate_distribution ')])
    with open(os.path.abspath(bppseqgenOptionsFile)) as f:
        bppseqgenModelLines = sorted([" ".join(x.split()) for x in f.readlines() if x.startswith('model ') or x.startswith('rate_distribution ')])
    assert coalhmmModelLines == bppseqgenModelLines, "model and rate distributions lines not identical in coalhmm and bppseqgen option files"

    import imp
    SimSpec = imp.new_module("SimSpec")
    exec open(inp).read() in SimSpec.__dict__
    sys.modules[SimSpec] = "SimSpec"

    def spec(**args):
        c1 = args["N1"] / args["NeRef"]
        c2 = args["N2"] / args["NeRef"]
        c3 = args["N3"] / args["NeRef"]
        c12 = args["N12"] / args["NeRef"]
        c123 = args["N123"] / args["NeRef"]
        t1 = args["T1"] / args["g"] / (2*args["NeRef"])
        t12 = args["T12"] / args["g"] / (2*args["NeRef"])
        t123 = args["T123"] / args["g"] / (2*args["NeRef"])
        speciesTree = P(c123, M(t12, [P(c3, S(1)), P(c12, M(t1, [P(c1, S(1)), P(c2, S(1))]))]))
        return speciesTree

    # simulate:
    simHook = CoaSimSimulationHook()
    if "recMap" in dir(SimSpec):
        seq = simulateCoasim(spec, length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u, addOutGroup=SimSpec.addOutGroup, \
                             T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123, \
                             N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123, \
                             optionsFile=os.path.abspath(bppseqgenOptionsFile), hook=simHook, recMap=SimSpec.recMap)
    else:
        seq = simulateCoasim(spec, length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u, addOutGroup=SimSpec.addOutGroup, \
                             T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123, \
                             N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123, \
                             optionsFile=os.path.abspath(bppseqgenOptionsFile), hook=simHook)

    assert len(simHook.recombinationPoints) == len(simHook.recombinationTimes), "%s vs %s" % (len(simHook.recombinationPoints), len(simHook.recombinationTimes))
    assert len(simHook.recombinationPoints) == len(simHook.trees) - 1, "%s vs %s" % (len(simHook.recombinationPoints), len(simHook.trees) - 1)


    # mapping between externa leaf label and state corresponding to the tree. note that
    # the simulated states do not distinguish state 0 and 1. so some 0s are really 1s
    stateMap = {"0": 0, "1": 3, "2": 2}

    # filter the recombination events to extract recombinations resulting in diffent topologies:
    recPoints = list()
    recTimes = list()
    recLeaves = list()
    fromState = list()
    toState = list()
    simStates = list()
    isSameTree = list()
    isSameTopology = list()
    isSameState = list()
    allRecPoints = list()
    allRecTimes = list()
    allRecLeaves = list()
    #simCoalTimes = list()

    def getTreeExternalLeaf(tree):
        left, right = tree.get_edges()
        if isinstance(left[0], Leaf):
            return left[0].identifier
        else:
            return right[0].identifier

    def getCoalTimes(tree):
        th = TreeHeight()
        tree.dfs_traverse(th)
        return sorted(list(set(th.max - x for x in th.depths)))


    for i in range(len(simHook.trees)-1):
        leftExtLeaf, rightExtLeaf = getTreeExternalLeaf(simHook.trees[i]),  getTreeExternalLeaf(simHook.trees[i+1])

        #######################################################
        leftCoalTimes = getCoalTimes(simHook.trees[i])
        rightCoalTimes = getCoalTimes(simHook.trees[i+1])
        #simCoalTimes.append([x * SimSpec.g * 2 * SimSpec.NeRef for x in sorted(leftCoalTimes)])
        t12 = SimSpec.T12 / SimSpec.g / (2*SimSpec.NeRef)
        leftState = stateMap[leftExtLeaf]
        if leftState == 0 and min(leftCoalTimes) > t12:
            leftState = 1
        rightState = stateMap[rightExtLeaf]
        if rightState == 0 and min(rightCoalTimes) > t12:
            rightState = 1
        if leftState != rightState:
            recPoints.append(simHook.recombinationPoints[i])
            recTimes.append(simHook.recombinationTimes[i])
            if coaSimLeafHack:
                # commented out because this functionality has a memory leak in CoaSim
                recLeaves.append(simHook.recombinationLeaves[i])
            fromState.append(leftState)
            toState.append(rightState)
        #######################################################


        isSameTree.append(int(str(simHook.trees[i]) == str(simHook.trees[i+1])))
        isSameTopology.append(int(leftExtLeaf == rightExtLeaf))
        isSameState.append(int(leftState == rightState))
        allRecPoints.append(simHook.recombinationPoints[i])
        allRecTimes.append(simHook.recombinationTimes[i])
        if coaSimLeafHack:
            # commented out because this functionality has a memory leak in CoaSim
            allRecLeaves.append(simHook.recombinationLeaves[i])


    # run coalhmm on simulated sequences:
    estimHook = ILS16estimationHook(minProb, maxSpan)
    t = estimate_ils16(seq, length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u,
                       addOutGroup=SimSpec.addOutGroup, T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123,
                       N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123,
                       optionsFile=os.path.abspath(coalhmmOptionsFile), hook=estimHook)

    ## estimHook = CTMCestimationHook()
    ## t = estimate_3sp_ctmc(seq, r=r, g=g, u=u, T1=T1, T2=T12, N1=N1, nstates=10, hook=estimHook)
    ## t = estimate_untitled_3_all(seq, T1 * u, T12 * u, 1.0/(N1 * u * 2*g), r/(u*g), nstates=10)
    ## t = estimate_untitled_3_all(seq, T1 * u, T12 * u, 1.0/(N1 * u * 2*g), r/(u*g), nstates=20, suffix= "_%s_%d" % (tag, iteration))
    ## os.system("cp %s sequence1.fasta" % seq)

    if not coaSimLeafHack:
        recLeaves = [-1] * len(recTimes)
        allRecLeaves = [-1] * len(recLeaves)

    if t is not None:
        stats = { 'ilsBases': estimHook.ilsBases,
                  'non_ilsBases': estimHook.non_ilsBases,
                  'simRecPoints': recPoints,
                  'simRecTimes': recTimes,
                  'simRecLeaves': recLeaves,
                  'simFromState': fromState,
                  'simToState': toState,
                  #'simCoalTimes': simCoalTimes,
                  'infRecPoints': [estimHook.startCoord[i] + (estimHook.endCoord[i] - estimHook.startCoord[i])/2.0 for i in range(len(estimHook.startCoord))],
                  'infStartCoords': estimHook.startCoord,
                  'infEndCoords': estimHook.endCoord,
                  'infFromState': estimHook.fromState,
                  'infToState': estimHook.toState,
                  'allSimRecPoints': allRecPoints,
                  'allSimRecTimes': allRecTimes,
                  'allSimRecLeaves': allRecLeaves,
                  'isSameTree': isSameTree,
                  'isSameTopology': isSameTopology,
                  'isSameState': isSameState,
                  'simParameters': dict([(k, str(v[0])) for k, v in t.data.items() if k.startswith('_')]),
                  'infParameters': dict([(k, v[0].strip()) for k, v in t.data.items() if not k.startswith('_')]),
                  'seqLength': SimSpec.L }
    #               'simParameters': dict([(k, v) for k, v in t.data.items() if k.startswith('_')]),
    #               'infParameters': dict([(k, v) for k, v in t.data.items() if not k.startswith('_')])}

        with open(outp, 'w') as f:
            pickle.dump(stats, f)

#         os.system("cp " + seq + "* .")
#         print seq
#         sys.exit()

    for f in glob.glob(seq + "*"):
        os.unlink(f)



def runSimulationsWithCoaSimAndILSCTMC(inp, outp, bppseqgenOptionsFile):
    """
    Simulate a number of 
    """

    assert os.path.exists(os.path.abspath(bppseqgenOptionsFile))

    import imp
    SimSpec = imp.new_module("SimSpec")
    exec open(inp).read() in SimSpec.__dict__
    sys.modules[SimSpec] = "SimSpec"
    
    def spec(**args):
        c1 = args["N1"] / args["NeRef"]
        c2 = args["N2"] / args["NeRef"]
        c3 = args["N3"] / args["NeRef"]
        c12 = args["N12"] / args["NeRef"]
        c123 = args["N123"] / args["NeRef"]
        t1 = args["T1"] / args["g"] / (2*args["NeRef"])
        t12 = args["T12"] / args["g"] / (2*args["NeRef"])
        t123 = args["T123"] / args["g"] / (2*args["NeRef"])
        speciesTree = P(c123, M(t12, [P(c3, S(1)), P(c12, M(t1, [P(c1, S(1)), P(c2, S(1))]))]))
        return speciesTree

    # simulate:
    simHook = CoaSimSimulationHook()
    if "recMap" in dir(SimSpec):
        seq = simulateCoasim(spec, length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u, addOutGroup=SimSpec.addOutGroup, \
                             T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123, \
                             N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123, \
                             optionsFile=os.path.abspath(bppseqgenOptionsFile), hook=simHook, recMap=SimSpec.recMap)
    else:
        seq = simulateCoasim(spec, length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u, addOutGroup=SimSpec.addOutGroup, \
                             T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123, \
                             N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123, \
                             optionsFile=os.path.abspath(bppseqgenOptionsFile), hook=simHook)

    assert len(simHook.recombinationPoints) == len(simHook.recombinationTimes), "%s vs %s" % (len(simHook.recombinationPoints), len(simHook.recombinationTimes))
    assert len(simHook.recombinationPoints) == len(simHook.trees) - 1, "%s vs %s" % (len(simHook.recombinationPoints), len(simHook.trees) - 1)


    # mapping between externa leaf label and state corresponding to the tree. note that
    # the simulated states do not distinguish state 0 and 1. so some 0s are really 1s
    stateMap = {"0": 0, "1": 3, "2": 2}

    # filter the recombination events to extract recombinations resulting in diffent topologies:
    recPoints = list()
    recTimes = list()
    recLeaves = list()
    fromState = list()
    toState = list()
    simStates = list()
    isSameTree = list()
    isSameTopology = list()
    isSameState = list()
    allRecPoints = list()
    allRecTimes = list()
    allRecLeaves = list()
    #simCoalTimes = list()
    
    def getTreeExternalLeaf(tree):
        left, right = tree.get_edges()
        if isinstance(left[0], Leaf):
            return left[0].identifier
        else:
            return right[0].identifier

    def getCoalTimes(tree):        
        th = TreeHeight()
        tree.dfs_traverse(th)
        return sorted(list(set(th.max - x for x in th.depths)))


    for i in range(len(simHook.trees)-1):
        leftExtLeaf, rightExtLeaf = getTreeExternalLeaf(simHook.trees[i]),  getTreeExternalLeaf(simHook.trees[i+1])

        #######################################################
        leftCoalTimes = getCoalTimes(simHook.trees[i])
        rightCoalTimes = getCoalTimes(simHook.trees[i+1])
        #simCoalTimes.append([x * SimSpec.g * 2 * SimSpec.NeRef for x in sorted(leftCoalTimes)])
        t12 = SimSpec.T12 / SimSpec.g / (2*SimSpec.NeRef)
        leftState = stateMap[leftExtLeaf]
        if leftState == 0 and min(leftCoalTimes) > t12:
            leftState = 1
        rightState = stateMap[rightExtLeaf]
        if rightState == 0 and min(rightCoalTimes) > t12:
            rightState = 1
        if leftState != rightState:
            recPoints.append(simHook.recombinationPoints[i])
            recTimes.append(simHook.recombinationTimes[i])
            # commented out because this functionality has a memory leak in CoaSim
            #recLeaves.append(simHook.recombinationLeaves[i])
            fromState.append(leftState)
            toState.append(rightState)
        #######################################################


        isSameTree.append(int(str(simHook.trees[i]) == str(simHook.trees[i+1])))
        isSameTopology.append(int(leftExtLeaf == rightExtLeaf))
        isSameState.append(int(leftState == rightState))
        allRecPoints.append(simHook.recombinationPoints[i])
        allRecTimes.append(simHook.recombinationTimes[i])
        # commented out because this functionality has a memory leak in CoaSim
        #allRecLeaves.append(simHook.recombinationLeaves[i])

    # run coalhmm on simulated sequences:	
    estimHook = CTMCilsHook(L=SimSpec.L, winsize=SimSpec.winsize)
    t = estimate_ctmcILS(seq, length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u,
                         addOutGroup=SimSpec.addOutGroup, T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123,
                         N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123,
                         hook=estimHook)

    if t is not None:

        from MultiPurpose import TimeSeries as TS
        observedTreeChanges = zip( *TS.summaryStats(allRecPoints, binIdx=0, binSize=SimSpec.winsize, stats=lambda buf: TS.poissonRateWithCI(buf, binSize=SimSpec.winsize)) )[2]
        observedTopologyChanges = zip( *TS.summaryStats(recPoints, binIdx=0, binSize=SimSpec.winsize, stats=lambda buf: TS.poissonRateWithCI(buf, binSize=SimSpec.winsize)) )[2]

        for tree, topol, exp in zip(observedTreeChanges, observedTopologyChanges, estimHook.expectedTransitions):
            print tree, topol, exp
        
        stats = { 'simRecPoints': recPoints,
                  'simRecTimes': recTimes,
                  'simRecLeaves': recLeaves,
                  'simFromState': fromState,
                  'simToState': toState,
                  'allSimRecPoints': allRecPoints,
                  'allSimRecTimes': allRecTimes,
                  'allSimRecLeaves': allRecLeaves,
                  #'simCoalTimes': simCoalTimes,
                  'isSameTree': isSameTree,
                  'isSameTopology': isSameTopology,
                  'isSameState': isSameState,
                  'observedTreeChanges': observedTreeChanges,
                  'observedTopologyChanges': observedTopologyChanges,
                  'expectedTransitions': estimHook.expectedTransitions,
                  'expectedTransitionsWinsize': SimSpec.winsize,                  
                  'simParameters': dict([(k, str(v[0])) for k, v in t.data.items() if k.startswith('_')]),
                  'infParameters': dict([(k, str(v[0])) for k, v in t.data.items() if not k.startswith('_')]),
                  'seqLength': SimSpec.L }


        with open(outp, 'w') as f:
            pickle.dump(stats, f)

#         os.system("cp " + seq + "* .")
#         print seq
#         sys.exit()

    for f in glob.glob(seq + "*"):
        os.unlink(f)


def evaluateSimulations(inp, outp):
    """
    
    """
    benchmarkStatsFileName = outp

    tableDict = dict()
    tableDictRandom = dict()

    totalLength = 0

    specId = -1
    
    randomShuffle = False

    ## for i in sorted(inp):
    for i in sorted(inp, key=lambda x: int(re.sub('[^\d]', '', x))): # input sorted numerially (only on numbers)

        try:
            iterationId =  re.search(r'itr(\d+)_', i).group(1)
        except:
            iterationId = 'NA'


        m = re.search(r'spec(\d+)\.', i)
        if m:
            assert int(m.group(1)) > int(specId)
            specId =  m.group(1)
        else:
            specId = 'NA'
        
        with open(i) as f:
            stats = pickle.load(f)

#         ##################################################################################################
#         # HACK to get length of simulated sequence because it was not retained in sim results
#         # at the time of coding. You can delete this for future reruns.
#         # HACK
#         print "HACK: skipping spec 91......"
#         if specId == '91':
#             # THIS IS THE ONE THAT DID NOT COMPLETE
#             specFileName =  re.search(r'\.(chrspec\d+)\.', i).group(1)
#             with open("full_chr_simulations/hapmapmap.%s.chr2.py" % specFileName) as f:
#                 stats["seqLength"] = int(float(re.search(r'L\s*=\s*(\S+)', f.read()).group(1)))            
#             totalLength += stats["seqLength"]
#             continue
#         ##################################################################################################

        assert "seqLength" in stats, stats.keys()

#         if "seqLength" not in stats:
#             stats["seqLength"] = 0 # this is just a hack so we can use this function for evaluating standard non-full chromosome simulations too.

        totalLength += stats["seqLength"]
                
        # hack: simRecTimes are not retained using MaCS...
        if "simRecTimes" not in stats:
            stats["simRecTimes"] = ['NA'] * len(stats["simRecPoints"])

        def stateMap(x):
            ## l = [0, 0, 2, 3] # simulations do not distinguish state 0 and 1
            l = [0, 1, 2, 3] # simulations do not distinguish state 0 and 1
            return l[x]


        def match(x, y):
            return abs(x[0] - y[0])
        #gap = 2000
        gap = 1000000

#         def match(x, y):
#             if not (x[1] == y[1] or x[2] == y[2]):
#                 return 10000
#             else:            
#                 return abs(x[0] - y[0])
#         gap = 1000


        def alignSeries(simRecPoints, infRecPoints):

            needle = Needle.Needle(simRecPoints)
            traceback = needle.align(infRecPoints, match, gap)

            for i, (simIdx, infIdx) in enumerate(traceback):

                if simIdx is not None and infIdx is not None:
                    aligned = 1
                else:
                    aligned = 0

                if simIdx is not None:
                    true = 1
                else:
                    true = 0

                if simIdx is not None:
                    simRecPoint, simFromState, simToState, simRecTime, simRecLeaf = simRecPoints[simIdx]
                else:
                    simRecPoint, simFromState, simToState, simRecTime, simRecLeaf = 'NA', 'NA', 'NA', 'NA', 'NA'

                if infIdx is not None:
                    infRecPoint, infFromState, infToState, infStartCoord, infEndCoord = infRecPoints[infIdx]
                else:
                    infRecPoint, infFromState, infToState, infStartCoord, infEndCoord = 'NA', 'NA', 'NA', 'NA', 'NA'

                if simIdx is not None and infIdx is not None:
                    distToClosestSim = infRecPoints[infIdx][0] - simRecPoints[simIdx][0]
                else:
                    prevSimPos, nextSimPos = None, None                    
                    j = i
                    while j >= 0:
                        j -= 1
                        if traceback[j][0] is not None:
                            prevSimPos = traceback[j][0]
                            break
                    j = i
                    while j < len(traceback)-1:
                        j += 1
                        if traceback[j][0] is not None:
                            nextSimPos = traceback[j][0]
                            break

                    if simIdx is not None:
                        l = list()
                        if prevSimPos is not None:
                            l.append(simRecPoints[simIdx][0] - prevSimPos)
                        if nextSimPos is not None:
                            l.append(simRecPoints[simIdx][0] - nextSimPos)
                        distToClosestSim = min(l, key=abs)

                    if infIdx is not None:
                        l = list()
                        if prevSimPos is not None:
                            l.append(infRecPoints[infIdx][0] - prevSimPos)
                        if nextSimPos is not None:
                            l.append(infRecPoints[infIdx][0] - nextSimPos)
                        distToClosestSim = min(l, key=abs)

#                 if simIdx is not None and infIdx is not None:
#                     print simRecPoint, infRecPoint simToState, infToState, infToState
#                         infRecPoint, infFromState, infToState, infStartCoord, infEndCoord = 'NA', 'NA', 'NA', 'NA', 'NA'
#                 infRecPoint, infFromState, infToState, infStartCoord, infEndCoord = 'NA', 'NA', 'NA', 'NA', 'NA'
                        
                for k, v in stats['simParameters'].items():
                    tableDict.setdefault(k, []).append(v)

                for k, v in stats['infParameters'].items():
                    tableDict.setdefault(k, []).append(v)


                if simRecPoint == 'NA':
                    cumulSimRecPoint = 'NA'
                else:
                    cumulSimRecPoint = simRecPoint + totalLength-stats['seqLength']

                if infRecPoint == 'NA':
                    cumulInfRecPoint = 'NA'
                else:
                    cumulInfRecPoint = infRecPoint + totalLength-stats['seqLength']


                if infStartCoord == 'NA':
                    cumulInfStartCoord = 'NA'
                else:
                    cumulInfStartCoord =  infStartCoord + totalLength-stats['seqLength']

                if infEndCoord == 'NA':
                    cumulInfEndCoord = 'NA'
                else:
                    cumulInfEndCoord =  infEndCoord + totalLength-stats['seqLength']

                tmp = { "iterationId": iterationId,
                        "specId": specId,
                        "aligned": aligned,
                        "true": true,
                        "simRecPoint": simRecPoint,
                        "simRecLeaf": ''.join(map(str, simRecLeaf)),
                        "cumulSimRecPoint": cumulSimRecPoint,
                        "infRecPoint": infRecPoint,
                        "cumulInfRecPoint": cumulInfRecPoint,
                        "simFromState": simFromState,
                        "infFromState": infFromState,
                        "simToState": simToState,
                        "infToState": infToState,
                        "distToClosestSim": distToClosestSim,
                        "simRecTime": simRecTime,
                        "infStartCoord": infStartCoord,
                        "cumulInfStartCoord": cumulInfStartCoord,
                        "infEndCoord": infEndCoord,
                        "cumulInfEndCoord": cumulInfEndCoord,
                        "seqLength": stats['seqLength'],
                        "seqStart": totalLength-stats['seqLength'],
                        "random": int(randomShuffle),
                        "seqEnd": totalLength }

                for k, v in tmp.items():
                    tableDict.setdefault(k, []).append(v)


# I SHOULD 

#        print len(stats['infToState']), sum(x==y for x,y in zip(map(stateMap, stats['infFromState']), map(stateMap, stats['infToState'])))


        simRecPoints = zip(stats['simRecPoints'], stats['simFromState'],
                           stats['simToState'], stats['simRecTimes'], stats['simRecLeaves'])
        # NOTE: here we delete 0 to 0 inferred events because these are not represented in simulated events we compare to
        infRecPoints = [t for t in zip(stats['infRecPoints'], map(stateMap, stats['infFromState']),
                                       map(stateMap, stats['infToState']), stats['infStartCoords'],
                                       stats['infEndCoords']) if t[1] != t[2]]

        if not len(infRecPoints):
            print "WARNING: DID NOT INFER ANY RECOMBINATIONS - SKIPPING"
            continue

        if not len(simRecPoints):
            print "WARNING: DID NOT SIMULATE ANY RECOMBINATIONS - SKIPPING"
            continue

        randomShuffle = False
        
        alignSeries(simRecPoints, infRecPoints)

        randomShuffle = True

        # get intervals between points
        infP, infF, infT, infS, infE = zip(*infRecPoints)
        intervals = [y-x for x, y in zip([0.0] + list(infP[0:-1]), infP)]
        # assign to states
        intervalsDict = dict()
        for i, t in zip(intervals, infF):
            intervalsDict.setdefault(t, []).append(i)
        #intervalsDict[infT[-1]].append(stats["seqLength"] - infP[-1])
        intervalsDict.setdefault(infT[-1], []).append(stats["seqLength"] - infP[-1])

        # random indexes
        indexes = range(len(infRecPoints))
        random.shuffle(indexes)
        # build new list
        newList = list()
        pos = 0
        for i in indexes:
            f, t = infF[i], infT[i]
            try:
                d = intervalsDict[f].pop(random.randint(0, len(intervalsDict[f])-1))
            except ValueError:
                # if len(intervalsDict[f])-1 is zero, representation may make it
                # diffferent from the other zero, whicch makes it die...
                if len(intervalsDict[f])-1 == 0:
                    d = intervalsDict[f].pop(0)
                else:
                    raise ValueError
            pos += d
            newList.append( (pos, f, t, 'NA', 'NA') )

        assert len(infRecPoints) == len(newList)
        infRecPoints = newList

        randomShuffle = True
        alignSeries(simRecPoints, infRecPoints)
        randomShuffle = False

        
    keys = sorted(tableDict.keys())
    values = [tableDict[k] for k in keys]
    with open(benchmarkStatsFileName, 'w') as f:
        print >>f, "\t".join(keys)
        for i in range(len(values[0])):
            print >>f, "\t".join(map(str, (values[j][i] for j in range(len(keys)))))


