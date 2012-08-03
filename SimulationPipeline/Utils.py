
import sys, os, glob
import CoaSim
from CoaSim.popStructure import Population as P, Sample as S, Merge as M, Migration as Mi, Growth
from SimulationPipeHotSpotCTMC import *
from newick.tree import Leaf
from CoalhmmPipeline import Table
import Needle
from optparse import OptionParser

def runSimulationsWithCoaSimScript():

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


def runSimulationsWithMaCS(inp, outp, minProb, maxSpan, coalhmmOptionsFile, ):
    """
    Simulate a number of 
    """

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
                           optionsFile=os.path.abspath("bppSeqGen.options"), hook=simHook, recMap=SimSpec.recMap)
    else:
        seq = simulateMaCS(length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u, addOutGroup=SimSpec.addOutGroup, \
                           T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123, \
                           N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123, \
                           optionsFile=os.path.abspath("bppSeqGen.options"), hook=simHook)

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
        leftExtLeaf, rightExtLeaf = getTreeExternalLeaf(simHook.trees[i]),  getTreeExternalLeaf(simHook.trees[i+1])
#         print 'left', simHook.trees[i], leftExtLeaf
#         print 'right', simHook.trees[i+1], rightExtLeaf
        if leftExtLeaf != rightExtLeaf:
            recPoints.append(simHook.recombinationPoints[i])
            fromState.append(stateMap[leftExtLeaf])
            toState.append(stateMap[rightExtLeaf])

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

    os.system("cp " + seq + "* tmp_macs/.")

    for f in glob.glob(seq + "*"):
        os.unlink(f)


def runSimulationsWithCoaSim(inp, outp, minProb, maxSpan, coalhmmOptionsFile):
    """
    Simulate a number of 
    """

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
                             optionsFile=os.path.abspath("bppSeqGen.options"), hook=simHook, recMap=SimSpec.recMap)
    else:
        seq = simulateCoasim(spec, length=SimSpec.L, NeRef=SimSpec.NeRef, r=SimSpec.r, g=SimSpec.g, u=SimSpec.u, addOutGroup=SimSpec.addOutGroup, \
                             T1=SimSpec.T1, T12=SimSpec.T12, T123=SimSpec.T123, \
                             N1=SimSpec.N1, N2=SimSpec.N2, N3=SimSpec.N3, N12=SimSpec.N12, N123=SimSpec.N123, \
                             optionsFile=os.path.abspath("bppSeqGen.options"), hook=simHook)

    assert len(simHook.recombinationPoints) == len(simHook.recombinationTimes), "%s vs %s" % (len(simHook.recombinationPoints), len(simHook.recombinationTimes))
    assert len(simHook.recombinationPoints) == len(simHook.trees) - 1, "%s vs %s" % (len(simHook.recombinationPoints), len(simHook.trees) - 1)


    # mapping between externa leaf label and state corresponding to the tree. note that
    # the simulated states do not distinguish state 0 and 1. so some 0s are really 1s
    stateMap = {"0": 0, "1": 3, "2": 2}

    # filter the recombination events to extract recombinations resulting in diffent topologies:
    recPoints = list()
    recTimes = list()
    fromState = list()
    toState = list()
    simStates = list()
    isSameTree = list()
    isSameTopology = list()
    allRecPoints = list()
    allRecTimes = list()

    def getTreeExternalLeaf(tree):
        left, right = tree.get_edges()
        if isinstance(left[0], Leaf):
            return left[0].identifier
        else:
            return right[0].identifier

    for i in range(len(simHook.trees)-1):
        leftExtLeaf, rightExtLeaf = getTreeExternalLeaf(simHook.trees[i]),  getTreeExternalLeaf(simHook.trees[i+1])

        isSameTree.append(int(str(simHook.trees[i]) == str(simHook.trees[i+1])))
        isSameTopology.append(int(leftExtLeaf == rightExtLeaf))
        allRecPoints.append(simHook.recombinationPoints[i])
        allRecTimes.append(simHook.recombinationTimes[i])

        if leftExtLeaf != rightExtLeaf:
            recPoints.append(simHook.recombinationPoints[i])
            recTimes.append(simHook.recombinationTimes[i])
            fromState.append(stateMap[leftExtLeaf])
            toState.append(stateMap[rightExtLeaf])

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

    if t is not None:
        stats = { 'simRecPoints': recPoints,
                  'simRecTimes': recTimes,
                  'simFromState': fromState,
                  'simToState': toState,
                  'infRecPoints': [estimHook.startCoord[i] + (estimHook.endCoord[i] - estimHook.startCoord[i])/2.0 for i in range(len(estimHook.startCoord))],
                  'infStartCoords': estimHook.startCoord,
                  'infEndCoords': estimHook.endCoord,
                  'infFromState': estimHook.fromState,
                  'infToState': estimHook.toState,
                  'allSimRecPoints': allRecPoints,
                  'allSimRecTimes': allRecTimes,
                  'isSameTree': isSameTree,
                  'isSameTopology': isSameTopology,
                  'simParameters': dict([(k, str(v[0])) for k, v in t.data.items() if k.startswith('_')]),
                  'infParameters': dict([(k, v[0].strip()) for k, v in t.data.items() if not k.startswith('_')]),
                  'seqLength': SimSpec.L }
    #               'simParameters': dict([(k, v) for k, v in t.data.items() if k.startswith('_')]),
    #               'infParameters': dict([(k, v) for k, v in t.data.items() if not k.startswith('_')])}

        with open(outp, 'w') as f:
            pickle.dump(stats, f)

        #os.system("cp " + seq + "* tmp_coasim/.")

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
            l = [0, 0, 2, 3] # simulations do not distinguish state 0 and 1
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
                    simRecPoint, simFromState, simToState, simRecTime = simRecPoints[simIdx]
                else:
                    simRecPoint, simFromState, simToState, simRecTime  = 'NA', 'NA', 'NA', 'NA'

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
                    while j < len(simRecPoints)-1:
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


        simRecPoints = zip(stats['simRecPoints'], stats['simFromState'], stats['simToState'], stats['simRecTimes'])
        # NOTE: here we delete 0 to 0 inferred events because these are not represented in simulated events we compare to
        infRecPoints = [t for t in zip(stats['infRecPoints'], map(stateMap, stats['infFromState']), map(stateMap, stats['infToState']), stats['infStartCoords'], stats['infEndCoords']) if t[1] != t[2]]

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
        intervalsDict[infT[-1]].append(stats["seqLength"] - infP[-1])
        # random indexes
        indexes = range(len(infRecPoints))
        random.shuffle(indexes)
        # build new list
        newList = list()
        pos = 0
        for i in indexes:
            f, t = infF[i], infT[i]
            d = intervalsDict[f].pop(random.randint(0, len(intervalsDict[f])-1))
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


