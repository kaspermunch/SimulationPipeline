import CoaSim, os, tempfile, Table, sys, subprocess
import pickle as pickle
from CoaSim.popStructure import Population as P, Sample as S, Merge as M
from .newick import *

from numpy import arange



# coalhmm_exe = "~tth/local/bin/coalhmm --noninteractive=yes"
# bppseqgen_exe = "~tth/local/bin/bppseqgen --noninteractive=yes"
# bppseqgen_dir = None

# coalhmm_exe = "/usr/local/bin/coalhmm --noninteractive=yes"
# bppseqgen_exe = "/usr/local/bin/bppseqgen --noninteractive=yes"
# bppseqgen_dir = None

coalhmm_exe = "./coalhmm --noninteractive=yes"
coalhmm_dir = "./scripts/Simulation/coalhmmOSX"
#coalhmm_dir = "./coalhmm-bpp-17-11-2011/package"

bppseqgen_exe = "./bppseqgen --noninteractive=yes" # specify exe in current dir
bppseqgen_dir = "./scripts/Simulation/bppseqgenOSX"  # current dir to change to before execution

class SimulationHook(object):

    def __init__(self):
        self.recombinationTimes = []
        self.recombinationPoints = []
        self.trees = []

    def run(self, arg, trees, args):

        # loop over nodes in the arg and collect recombination points and times:
        recombinationEvents = dict()
        for n in arg.nodes:            
            if isinstance(n, CoaSim.RecombinationNode) and n.isAncestral(n.recombinationPoint):
                recombinationEvents[n.recombinationPoint] = n.eventTime

        # sort and scale points and times:
        for p, t in sorted(recombinationEvents.items()):
            self.recombinationPoints.append(int(round(p  * args["length"])))
            self.recombinationTimes.append(int(round(t * (2.0*args["NeRef"]) * args["g"])))

        # get the tree for each interval
        for interval in arg.intervals:
            self.trees.append(parse_tree(str(interval.tree)))

        # if we use recombination rate variation we need to use the adjusted interval
        # starts instead of the recombinationPoints. We rely on that trees are sorted too.
        if "recMap" in args:            
            treeStartList = list()
            for i, (start, end, tree) in enumerate(trees):
                if i != 0:
                    treeStartList.append(int(round(start  * args["length"])))
            self.recombinationPoints = treeStartList

#             f = open("rec.tbl", 'w')
#             for j in self.recombinationPoints:
#                 print >>f, j
#             f.close()
# 
#             f = open("tre.tbl", 'w')
#             for j in treeStartList:
#                 print >>f, j
#             f.close()
            
#             arg.intervals[i].start = start
#             arg.intervals[i].end = end
# 
#             if i != 0:
#                 recPoints.append(start)
#         recPoints.sort()
#         origIdx = [t[1] for t in sorted(zip([n.recombinationPoints for n in arg.nodes], range(len(arg.nodes))))]
#         for k, o in enumerate(origIdx):
#             if isinstance(CoaSim.RecombinationNode, arg.nodes[o]):
#                 arg.nodes[o].recombinationPoint = recPoints[k]


class ILSestimationHook(object):

    def __init__(self):
        self.startCoord = []
        self.endCoord = []
        self.fromState = list()
        self.toState = list()
        self.infStates = list()

    def run(self, prefix, args):
        prevInferredState = None
        prevState = None
        prevPos = None
        f = open(prefix + ".posterior.tbl", 'r')
        next(f) # remove header
        prevPos = 0
        pos = 0
        for l in f:
            lst = list(map(float, l.split()[1:5]))
            #argmax = lambda lst: max(izip(lst, xrange(len(lst))))[1]
            maxProb = max(lst)
            self.infStates.append(lst.index(maxProb))
            if maxProb > 0.5:
                state = lst.index(maxProb)
                if prevState is not None and state != prevState:
                    self.startCoord.append(prevPos)
                    self.endCoord.append(pos)
                    self.fromState.append(prevState)
                    self.toState.append(state)
                prevState, prevPos = state, pos
            pos += 1


class CTMCestimationHook(object):

    def __init__(self):
        self.startCoord = []
        self.endCoord = []
        self.fromState = list()
        self.toState = list()
        self.infStates = list()

    def run(self, prefix, args):
        f = open("%s_model.pickle")
        model = pickle.load(f)
        f.close()
        tree_map = model.treemap
        # PARSE MODEL TO SO WE CAN DISTINGUISH MODELS WITH DIFFERENT TOPOLOGIES AND WIGHIN THOSE WHAT THE BRANCHLENGTHS (IN TERMS OF BINS) ARE....

        prevInferredState = None
        prevState = None
        prevPos = None
        f = open(prefix + "_pds.txt", 'r')
        next(f) # remove header
        prevPos = 0
        pos = 0
        for l in f:
            lst = list(map(float, l.split()[1:5]))
            #argmax = lambda lst: max(izip(lst, xrange(len(lst))))[1]
            maxProb = max(lst)
            self.infStates.append(lst.index(maxProb))
            if maxProb > 2.0 / args["nstates"]:
                state = lst.index(maxProb)
                if prevState is not None and state != prevState:
                    self.startCoord.append(prevPos)
                    self.endCoord.append(pos)
                    self.fromState.append(prevState)
                    self.toState.append(state)
                prevState, prevPos = state, pos
            pos += 1

class TreeHeight(tree.TreeVisitor):
    def __init__(self):
        self.height = 0
        self.max = 0

    def reset(self):
        TreeHeight.__init__(self)
    def result(self):
        return self.max
    def get_result(self):
        return self.max

    def pre_visit_edge(self,src,bootstrap,length,dst):
        self.height = self.height + length
        self.max = max(self.max, self.height)
        
    def post_visit_edge(self,src,bootstrap,length,dst):
        self.height = self.height - length

class LeafCount(tree.TreeVisitor):
    def __init__(self):
        self.count = 0

    def visit_leaf(self,l):
        self.count = self.count + 1
        

class TreeScale(tree.TreeVisitor):
    def __init__(self, factor):
        self.factor = factor
        self.parent = []
        self.top = None

    def pre_visit_tree(self,t):
        self.parent.append(tree.Tree())
        
    def post_visit_tree(self,t):
        self.top = self.parent.pop()
        
        
    def post_visit_edge(self,src,bootstrap,length,dst):
        t = self.parent.pop()
        self.parent.append(t)
        t.add_edge((self.top, bootstrap, length*self.factor)) 

    def visit_leaf(self,l):
        self.top = tree.Leaf(l.identifier)

def get_tau(n, **args):
    name = "T"
    for i in range(1, n):
        name = name + str(i)

    if name not in args:
        raise Exception("Parameter '" + name + "' must be specified")

    T = args[name]

    return T / args["g"] / (2 * args["NeRef"])

def addOutGroup(trees, **args):
    oldtrees = trees
    trees = []
    th = TreeHeight()
    lc = LeafCount()
    
    for (s,e,t) in oldtrees:
        th.max = 0
        lc.count = 0
        t.dfs_traverse(th)
        t.dfs_traverse(lc)

        new_height = get_tau(lc.count +1, **args)

        l = tree.Leaf(str(lc.count))
        t2 = tree.Tree()
        t2.add_edge((t, 0, new_height - th.max))
        t2.add_edge((l, 0,new_height))
        trees.append((s,e,t2))

   

    return trees

def scale_trees(trees, factor):
    ts = TreeScale(factor)
    oldtrees = trees
    trees = []
    for (s,e,t) in oldtrees:
        t.dfs_traverse(ts)
        trees.append((s,e,ts.top))

    return trees

def bppSeqGen(trees, **args):
    #write tree to tmp file
    (fd1, path1) = tempfile.mkstemp();
    f1 = os.fdopen(fd1, "w")
    for (s,e,t) in trees:
        f1.write(str(s) + " " + str(e) +" " + str(t) + ";\n")

#    print path1

    f1.close()
    
    #run bppseqgen
    (fd2, path2) = tempfile.mkstemp();
    f2 = os.fdopen(fd2, "w")

    logfilestd = ""
    if "logFile" in args:
        logfilestd = " > " + args["logFile"]
    
    cmd = bppseqgen_exe + " param=" + str(args["optionsFile"]) + " number_of_sites=" + \
        str(args["length"]) + " input.tree.file=" + path1 + " output.sequence.file=" + path2 + \
              logfilestd

#     print cmd
#    os.system(cmd)

    p = subprocess.Popen(cmd, env=os.environ, shell=True, cwd=bppseqgen_dir, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    stdout, stderr = p.communicate()
#     print stdout
#     print stderr

    #return path to sequence
    #todo: remove tmp file
    return path2


def normalizeSegments(segments):
#    print segments
    start, end = segments[0][0], segments[-1][1]
    length = end-start
    t = segments[0]
    norm_segments = [(0, (t[1]-start)/float(length), t[2])]
    for t in segments[1:]:
        norm_segments.append((((t[0]-start)/float(length)), (t[1]-start)/float(length), t[2]))
#    print norm_segments
    return norm_segments


def generateMapping(rRateSegments):
    """
    Generates the mapping from a specification of recombination rates on segmnets from 0 to 1.
    """
    segments = normalizeSegments(rRateSegments)
    mapList = list()
    meanRate = float(sum([x[2] * (x[1] - x[0]) for x in segments]))
    meanRateRest = meanRate
    totLength = 0.0
    while len(segments) > 1:
        start, end, rate = segments.pop(0)
        meanRateRest = (meanRate - (end-start)*rate)/(1-(end-start))
        oldLength = end-start
        if rate == 0:
            newLength = 0.0
            f = float('inf')
        else:
            newLength = 1/((((1-oldLength)/oldLength) * (meanRateRest/rate)) + 1)
            f = oldLength / newLength
        totLength += newLength
#        mapList.append((totLength, f))
        mapList.append((totLength, oldLength, newLength))
    start, end, rate = segments.pop(0)
    meanRateRest = (meanRate - (end-start)*rate)/(1-(end-start))
    oldLength = end-start
    if rate == 0:
        newLength = 0.0
        f = float('inf')
    else:
        newLength = 1/((((1-oldLength)/oldLength) * (meanRateRest/rate)) + 1)
        f = oldLength / newLength
    totLength += newLength
#     mapList.append((totLength, f))
    mapList.append((totLength, oldLength, newLength))

    # there is some numerical impression somewhere that means that the last coordinate is
    # not exactly 1 if the map is very large. so we just make sure it is very close, and
    # then set it to one:
    assert abs(1 - mapList[-1][0]) < 1e-10
#     mapList[-1] = (1.0, mapList[-1][1])
    mapList[-1] = (1.0, mapList[-1][1], mapList[-1][2])

    return mapList


def mapIntervalsAndParseTrees(intervals, segments):
    """
    Maps from simulated intervals to those resulting from varying recombination rate as
    specified by mapList.
    """
    mapList = generateMapping(segments)
    j = 0    
    prevEnd = 0
    newIntervals = list()
    for i, interval in enumerate(intervals):
        interval_start, interval_end = interval.start, interval.end
        end = prevEnd
        length = interval_end - interval_start
        while True:
            mapStart = 0
            if j > 0:
                mapStart = mapList[j-1][0]
            if interval_start + length > mapList[j][0]:
                if mapList[j][2] == 0:
                    end += mapList[j][1]
                else:
                    end += (mapList[j][0] - max(interval_start, mapStart)) * (mapList[j][1] / mapList[j][2])
                j += 1
            else:
                if mapList[j][2] == 0:
                    end += mapList[j][1]
                else:
                    end += (interval_end - max(interval_start, mapStart)) * (mapList[j][1] / mapList[j][2])
                break

        start = prevEnd
        prevEnd = end
        newIntervals.append((start, end, parse_tree(str(interval.tree))))

    # this handles the special case where the last specified recombination rate segment
    # has rate zero:
    if mapList[-1][2] == 0:
        newIntervals[-1] = (newIntervals[-1][0], newIntervals[-1][1] + mapList[-1][1], newIntervals[-1][2])

#     with open("test.tbl", 'w') as f:
#         for start, end, tree in newIntervals:
#             print >>f, start + (end-start)/2.0
#     sys.exit()
    
    return newIntervals

# 
# # List of tuples of start, end, rate of segments from 0 to 1. Rates are only relative so
# # mean rate must must correspond to the one used for simulation.
# rRateSegments = [(1000, 1250, 20), (1250, 2000, 10)]
# #rRateSegments = [(0, 0.25, 20), (0.25, 1, 10)]
# #rRateSegments = [(0, 0.75/2, 10), (0.75/2, 1-(0.75/2), 20), (1-(0.75/2), 1, 10)]
#     
# 
# 
# # Dummy intervals for testing
# by = 0.1
# starts = arange(0,1+by,by)
# intervals = zip(starts[0:-1], starts[1:])
# 
# 
# # runs like this
# normSegments = normalizeSegments(rRateSegments)
# mapList = generateMapping(normSegments)
# newIntervals = mapIntervals(intervals, mapList)
# 
# #
# # newIntervals = mapIntervals(intervals, rRateSegments)
# 
# 
# print newIntervals

def simulateCoasim(spec, **args):
    #test if all args are there
    if "length" not in args:
        raise Exception("Sequence 'length' parameter required")
    
    if "NeRef" not in args:
        raise Exception("'NeRef' parameter required")
    
    if "r" not in args:
        raise Exception("'r' (recombination rate?) parameter required")

    if "g" not in args:
        raise Exception("'g' (Generation time) parameter required")

    if "u" not in args:
        raise Exception("'u' (mutation rate) parameter required")

    if "optionsFile" not in args:
        raise Exception("'optionsFile' parameter required")

    if "migration" not in args:
        args["migration"] = []

    #test T's and N's
    count = 0;
    t = "T1"
    while t in args:
        count = count +1
        t = t + str(len(t))

    #Test N
    n = "N"
    for i in range(1, count+1):
        if "N" + str(i) not in args:
            raise Exception("'N" + str(i) + "' parameter required")
        n = n + str(i)
        if n not in args:
            raise Exception("'" + n + "' parameter required")

     
    #run the suff
    
    st = spec(**args)
    kei = True
    if "keepEmptyIntervals" in args:
        kei = args["keepEmptyIntervals"];
    
#    print "simulating"
        
    arg = CoaSim.simulate([], st, rho = 4 * args["NeRef"] * args["length"] * args["r"],
            keepEmptyIntervals=kei, migrationSpec=args["migration"]) #what about 4
    
    if "recMap" in args:
        # map intervals and parse trees
        trees = mapIntervalsAndParseTrees(arg.intervals, args["recMap"])
    else:
        # parse trees
        trees = []
        for interval in arg.intervals:
            trees.append((interval.start, interval.end, parse_tree(str(interval.tree))))        

    #reroot trees
    if "addOutGroup" in args and args["addOutGroup"]:
        trees = addOutGroup(trees, **args)
    
    #scale trees
    trees = scale_trees(trees, 2.0*args["NeRef"]*args["g"]*args["u"])

    if "tree_height_file" in args:
        fout = open(args["tree_height_file"], "w")
        bps = args["length"]
        if "tree_visitor" not in args:
            th = TreeHeight()
        else:
            th = args["tree_visitor"]
        for (start, end, tree) in trees:
            th.reset();
            tree.dfs_traverse(th)
            h = th.get_result()
            to_write = "%s\n" % str(h)
            for i in range(int(start*bps), int(end*bps)):
                fout.write(to_write)
        fout.close()

    #for (s,e,t) in trees:
        #print s,e,t

    if "hook" in args:
        args["hook"].run(arg, trees, args)
#         args["hook"].parseTrees(arg, args)
    
    return bppSeqGen(trees, **args)


def rerootMacsTree(tree):
    s = str(tree)
    s.replace(')', ' ')
    s = s.replace(')', ' ')
    s = s.replace('(', ' ')
    s = s.replace("'0'", ' ')
    s = s.replace("'1'", ' ')
    s = s.replace("'2'", ' ')
    s = s.replace(':', ' ')
    s = s.replace(';', ' ')
    s = s.replace(',', ' ')
    s = s.replace('  ', ' ')
    s = s.replace('  ', ' ')
    s = s.replace('  ', ' ')
    s = s.strip()
    ss = s.split()
    max = 0
    for e in ss:
        if float(e) > max:
            max = float(e)
        
    s = str(tree)
    s = s.replace(';', '')
    val = tau123 - max
    s = '(' + s + ' : ' + str(val) + " , '3' : " + str(tau123) + ' );'
    return s


def simulateMacs(**args):

    outputdir = tempfile.mkdtemp()

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
        
    cmd = macs_exe + " 3 " + str(args["L"]) + " -T -r " + str(4 * args["NeRef"] * args["r"]) + " -I 3 1 1 1 " + \
    " -n 1 " + str(args["N1"] / args["NeRef"]) + " -n 2 " + str(args["N1"] / args["NeRef"]) + " -n 3 " + str(args["N1"] / args["NeRef"]) +  \
    " -ej " + str(args["tau1"]) + " 2 3 -en " + str(args["tau1"]) + " 2 " + str(args["N12"] / args["NeRef"]) +  " -en " + str(args["tau1"]) + " 3 " + str(args["N12"] / args["NeRef"]) + \
    " -ej " + str(args["tau12"] ) +" 3 1 -en " + str(args["tau12"]) + " 1 " + str(args["N123"] / args["NeRef"]) 
    print(cmd)
    os.system(cmd + " > " + outputdir + "/macs.out 2> " + outputdir +"/macs.err")
    
    trees = list()
    
    treelengthmatcher = re.compile(': *-?\d*\.?\d*e?-?\d*')
    macsmatch = re.compile('.*\[(\d*)\](.*)')

    linecount = 0;
    out = open(outputdir + "/trees.dnd", "w") 
    fin = open(outputdir + "/macs.out", "r")
    tmptree = open("tmptree" , "w")
    last = 0;
    for item in fin :
        linecount = linecount +1
        #first 2 lines does not contain trees
        if(linecount < 3):
            continue;
        
        m = macsmatch.match(item)
        l = int(m.group(1))
        tree = m.group(2)
        tree = tree.replace("0:", "'0':")
        tree = tree.replace("1:", "'1':")
        tree = tree.replace("2:", "'2':")
        tree = rerootMacsTree(tree)
        
        start = float(last)/L
        end = float(last+l)/L

        lens = treelengthmatcher.findall(tree)
        for le in lens :
            fl = float(le.replace(':', '')) * 2*NeRef*u*g
            tree = tree.replace(le, ':' + str(fl))

        trees.append((start, end, parse_tree(tree)))
        
        last = last + l

    outputdir

    return bppSeqGen(trees)




def estimate_ils09(sequence, **args):
    #test if all args are there
    if "NeRef" not in args:
        raise Exception("'NeRef' parameter required")
    
    if "r" not in args:
        raise Exception("'r' (recombination rate?) parameter required")

    if "g" not in args:
        raise Exception("'g' (Generation time) parameter required")

    if "u" not in args:
        raise Exception("'u' (mutation rate) parameter required")

    if "optionsFile" not in args:
        raise Exception("'optionsFile' parameter required")

    #test T's and N's
    count = 0;
    t = "T1"
    while t in args:
        count
        t = t + str(len(t))

    #Test N
    n = "N"
    for i in range(1, count+1):
        if "N" + str(i) not in args:
            raise Exception("'N" + str(i) + "' parameter required")

    
    cmd = coalhmm_exe + " param=" + args["optionsFile"]
    
    NeRef = args["NeRef"]
    g = args["g"]
    u = args["u"]
    r = args["r"]
    T1 = args["T1"] / g / (2*NeRef)
    T12 = args["T12"] / g / (2*NeRef)
    T123 = args["T123"] / g / (2*NeRef)
    N1 = args["N1"] 
    N12 = args["N12"]
    N123 = args["N123"]

    param = "tau1=" + str(T1 * u * 2 * NeRef * g) + " tau2=" + str((T12-T1) * u * 2 * NeRef * g) + " c2=" +\
        str(((T123- T12) * NeRef - 4/3 * N123) * u * 2 * g) + " theta1=" + str(N12 * u * 2 * g) + " theta2=" + \
        str(N123 * u * 2 * g) + " rho=" + str(r / (u * g)) + " theta12=" + str(N12 * u * 2 * g) + " DATA=" + sequence

#     p = subprocess.Popen(cmd + " " + param, shell=True, cwd=coalhmm_dir)
#     p.wait()
#     print cmd + " " + param
    p = subprocess.Popen(cmd + " " + param, env=os.environ, shell=True, cwd=coalhmm_dir, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    stdout, stderr = p.communicate()
    print(stdout)
    print(stderr)
    if "Oups... something abnormal happened!" in stderr or "Optimization failed because likelihood function returned NaN" in stderr:
        print("estimation failed")
        return None


    if "hook" in args:
        args["hook"].run(sequence, args)

    input_t = Table.Table().add_row(_tau1=T1 * u * 2 * NeRef * g, _tau2=(T12-T1) * u * 2 * NeRef * g, \
                              _c2=((T123- T12) * NeRef - 4/3 * N123) * u * 2 * g, _theta1=N12 * u * 2 * g, \
                              _theta2=N123 * u * 2 * g, _rho=r / (u * g), _theta12=N12 * u * 2 * g)

    return Table.Table().load_kv(sequence + ".user.txt").join(input_t)
 
# # def estimate_untitled_3_all(seq, T1, T2, C, R, nstates=5):
# #     def logL_all(model, obs, colmap, c):
# #         if min([c]) <= 0:
# #             return -1e18
# # 	print "C =", c
# #         return logLikelihood(model, obs, colmap, [c,c], [R,R], [0.0,0.0], [0.0,T1])
# # 
# #     #antal epoker i model 2, hvad de konvergerer til 0, antal states i de epoker hhv. 1, nstates
# #     ssplit = build_epoch_seperated_model(2, [[0,0]], [1,nstates], None)
# #     names = ["'0'", "'1'"]
# #     est = optimize_f(seq, ssplit, names, logL_all, (C,))
# #     
# #     obs, colmap = readObservations(seq, names)
# #     _, pd, tbps, pi = logLikelihood(ssplit, obs, colmap, [est[0],est[0]], [R,R], [0.0,0.0], [0.0,T1], posterior_decoding=True)
# # 
# #     print pd.get_no_rows(), "x", pd.get_no_columns()
# # 
# #     f_bps = open("bps.txt", "w")
# #     print >>f_bps, "\t".join(map(str, tbps))
# #     f_bps.close()
# #     f_pi = open("pi.txt", "w")
# #     print >>f_pi, "\t".join([str(pi[x]) for x in xrange(len(pi))])
# #     f_pi.close()
# # 
# #     pd_file = open("pds.txt", "w")
# #     print >>pd_file, matrix_to_string(pd)
# #     pd_file.close()
# # 
# #     return Table.Table().add_row(C=est[0],R=R,tau_1=T1)
# 
# 
# 
# def print_matrix_to_file(M, f):
#     for row in xrange(M.get_no_rows()):
#         line = []
#         for col in xrange(M.get_no_columns()):
#             line.append(str(M[row,col]))
#         print >>f, ("\t".join(line))
# 
# def matrix_to_string(M):
#     res = []
#     for row in xrange(M.get_no_rows()):
#         sub_res = []
#         for col in xrange(M.get_no_columns()):
#             sub_res.append(str(M[row,col]))
#         res.append("\t".join(sub_res))
#     return "\n".join(res)
# 
# import sys
# sys.path.append( "/users/aeh/birc/working/" )
# sys.path.append( "/users/aeh/parredhmm/boostpython/" )
# sys.path.append( "/users/asand/public_html/" )
# #/users/aeh/birc/CoaSim/coasim-python-1.2/Python/build/lib.linux-i686-2.4/
# sys.path.append( "/users/aeh/birc/nymodel")
# import os
# 
# from model import build_epoch_seperated_model, build_simple_model
# from optimize import logLikelihood, readObservations
# from scipy.optimize import fmin
# 
# def optimize_f(model, obs, f_logL, init):
#     return fmin(lambda x: -f_logL(model, obs, *x), init)
# 
# def estimate_untitled_all_Nblocks(model, obs, T1, T2, C, R, outfile="/dev/null"):
#     def logL_all(model, all_obs, t1, t2, c, r):
#         if min([t1,t2,c,r]) <= 0 or t2 <= t1:
#             return -1e18
#         res = 0
#         for obs, colmap in all_obs:
#             res += logLikelihood(model, obs, colmap, [c]*3, [r]*3, [0.0]*3, [0.0,t1,t2])
#         os.system("echo '%s' >> %s" % ('\t'.join(map(str, [t1, t2, c, r, res])), outfile))
#         return res
# 
#     os.system("echo '%s' > %s" % ('\t'.join(map(str, ["t1", "t2", "c", "r", "logL"])), outfile))
#     est = optimize_f(model, obs, logL_all, (T1,T2,C,R))
# 
#     #       t1,     t2,     c,      r
#     return (est[0], est[1], est[2], est[3])
# 
# def estimate_3sp_ctmc(sequence, **args):
# 
# #     NeRef = N1 = N2 = N3 = 30000.0
# #     N12 = 30000.0
# #     N123 = 30000.0
# #     g = 20
# #     T1 = 0.5e6 
# #     T12 = 5e6
# #     L = 500e3 # 500 kbp
# #     r = 1.5e-8
# #     u = 1e-9
#     N1 = args["N1"]
#     T1 = args["T1"]
#     T2 = args["T2"]
#     r = args["r"]
#     u = args["u"]
#     g = args["g"]
#     nstates = args["nstates"]
# 
#     prefix = "ctmc_test"
# 
# #     if args["N1"] != args["N2"] or args["N1"] != args["N12"] or args["N1"] != args["N123"]:
# #         raise Exception("only one pop size allowed")
# 
#     model = build_epoch_seperated_model(3, [[0,0,1], [0,0]], [1,nstates,nstates], None)
# 
#     names = ["'0'", "'1'", "'2'"]
#     seqs = [sequence]
#     all_obs = []
#     for seq in seqs:
#         obs, colmap = readObservations(seq, names)
#         all_obs.append((obs, colmap))
# 
#     e_t1, e_t2, e_c, e_r = estimate_untitled_all_Nblocks(model, all_obs, T1 * u, T2 * u, 1.0/(N1 * u * 2*g), r/(u*g), outfile="%s_progress.txt" % prefix)
# 
#     obs, colmap = all_obs[0]
#     logL, pd, tbps, pi, pi_count, T_count, E_count = logLikelihood(model, obs, colmap, [e_c,e_c,e_c], [e_r,e_r,e_r], [0.0,0.0,0.0], [0.0,e_t1,e_t2], posterior_decoding=True)
# 
#     import cPickle as pickle
#     f = open("%s_model.pickle" % prefix, 'w')
#     pickle.dump(model, f)
#     f.close()
# 
#     f_bps = open("%s_bps.tbl" % prefix, "w")
#     print >>f_bps, "\t".join(map(str, tbps))
#     f_bps.close()
#     
#     f_pi = open("%s_pi.tbl" % prefix, "w")
#     print >>f_pi, "\t".join([str(pi[x]) for x in xrange(len(pi))])
#     f_pi.close()
# 
#     pd_file = open("%s_pds.txt" % prefix, "w")
#     print_matrix_to_file(pd, pd_file)
#     pd_file.close()
# 
#     if "hook" in args:
#         args["hook"].run(prefix, args)
# 
#     return(Table.Table().add_row(C=e_c,R=e_r,tau_1=e_t1,tau_2=e_t2))
