
#coaSimLeafHack = True
coaSimLeafHack = False

import sys
if coaSimLeafHack:
    sys.path.insert(0, '/Users/kasper/Desktop/coasim_trunk/Python/build/lib.macosx-10.5-x86_64-2.7')
import CoaSim, os, tempfile, subprocess, re, random
import cPickle as pickle
from CoaSim.popStructure import Population as P, Sample as S, Merge as M
from newick import *
from itertools import izip
from time import sleep
from numpy import arange

from CoalhmmPipeline import Table

from newick.tree import Leaf


# coalhmm_exe = "~tth/local/bin/coalhmm --noninteractive=yes"
# bppseqgen_exe = "~tth/local/bin/bppseqgen --noninteractive=yes"
# bppseqgen_dir = None

# coalhmm_exe = "/usr/local/bin/coalhmm --noninteractive=yes"
# bppseqgen_exe = "/usr/local/bin/bppseqgen --noninteractive=yes"
# bppseqgen_dir = None

coalhmm_exe = "./coalhmm --noninteractive=yes"
#coalhmm_dir = os.path.join(os.path.dirname(__file__), 'coalhmmOSX')
coalhmm_dir = "."

macs_exe = "./macs-patched"
macs_dir = "./scripts/Simulation/macsOSX"

bppseqgen_exe = "./bppseqgen --noninteractive=yes" # specify exe in current dir
#bppseqgen_dir = os.path.join(os.path.dirname(__file__), 'bppseqgenOSX')
bppseqgen_dir = "."


#
## cmtc hook:
#from coalhmm.optimize import default_bps
#from coalhmm.mini_hmm import calc_forward_backward
#
#from scipy import zeros, array, arange, int16, around
#from scipy.weave import inline
#from scipy import *

## run ctmc
#from runCTMC import runILSctmc

class CTMCilsHook(object):

    def __init__(self, **kwargs):
        vars(self).update(kwargs)
        self.expectedTransitions = list()

    def run(self, modelILS, COL_MAP, all_obs, est_c, est_r, est_m, est_t):

        def fill_table(F, T, E, n, m, S, obs):
            assert 0 <= n < m < len(obs)
            K = m - n
            assert S > 0 and K > 0
            assert F.shape == (S,)
            assert T.shape == (S, S)
            assert E.shape[0] == S

            Tnm_prev = zeros((S,K))
            Tnm = zeros((S,K))
            Tnm_prev[:,0] = F[:]

            for xm in obs[n:m]:
                Tnm_prev = Tnm_prev / Tnm_prev.sum()
                #print around(Tnm_prev, 4)
                #Tnm_prev = (Tnm_prev.T / Tnm_prev.sum(axis=1)).T
                for s in xrange(S):
                    for k in xrange(K):
                        x = Tnm_prev[s,k] * T[s,s]
                        if k > 0:
                            for sm in xrange(S):
                                if sm != s:
                                    x += Tnm_prev[sm, k-1]*T[sm,s]
                        Tnm[s,k] = x*E[s, xm]
                Tnm_prev,Tnm = Tnm,Tnm_prev
            Tnm_prev = Tnm_prev / Tnm_prev.sum()
            #print around(Tnm_prev, 4)
            #Tnm_prev = (Tnm_prev.T / Tnm_prev.sum(axis=1)).T
            return Tnm_prev

        def fill_table_fast(F, T, E, n, m, S, obs):
            assert 0 <= n < m < len(obs)
            K = m - n
            assert S > 0 and K > 0
            assert F.shape == (S,)
            assert T.shape == (S, S)
            assert E.shape[0] == S

            Tnm_prev = zeros((S,K))
            Tnm = zeros((S,K))
            Tnm_prev[:,0] = F[:]

            support_code = """
            void normalize_slice_inplace(double *slice, int S, int K) {
                int i;
                double sum = 0;
                for (i = 0; i < S*K; i++)
                    sum += slice[i];
                double norm = 1.0 / sum;
                for (i = 0; i < S*K; i++)
                    slice[i] *= norm;
            }
            void print_slice(double *slice, int S, int K) {
                int s,k;
                for (s = 0; s < S; s++)
                {
                    for (k = 0; k < K; k++)
                    {
                        printf("%0.4f\\t", slice[s*K + k]);
                    }
                    printf("\\n");
                }
            }
            """
            code = """
            #line 119 "rekomb.py"
            short *xmi = obs + n; // including n
            short *xme = obs + m; // excluding m
            for ( ; xmi != xme; xmi++ )
            {
                int xm = *xmi;
                normalize_slice_inplace(Tnm_prev, S, K);
                //print_slice(Tnm_prev, S, K);
                #pragma omp parallel for schedule(dynamic, 100)
                for (int k = 0; k < K; k++)
                {
                    for (int s = 0; s < S; s += 1)
                    {
                        double x = Tnm_prev[s*K + k] * T[s*S + s];
                        for (int sm = 0; k > 0 && sm < S; sm++)
                        {
                            if (sm != s)
                                x += Tnm_prev[sm*K + k-1] * T[sm*S + s];
                        }
                        Tnm[s*K + k] = x*E2(s, xm);
                    }
                }
                double *tmp = Tnm; Tnm = Tnm_prev; Tnm_prev = tmp;
            }
            normalize_slice_inplace(Tnm_prev, S, K);
            //print_slice(Tnm_prev, S, K);
            """

            inline(code, support_code=support_code,
                   arg_names=['Tnm','Tnm_prev','S','K','m','n','obs','T','E'],
                   ## extra_compile_args=["-fopenmp"],
                   ## extra_objects=["-fopenmp"]
                   )

            # Just doing pointer swaps in the C code won't change what the python
            # objects point to. So if there has been an odd number of swaps the
            # result is stored in Tnm, and otherwise in Tnm_prev
            if K % 2 == 1:
                return Tnm
            else:
                return Tnm_prev

        # matrices for forward backrward
        def get_matrices(model, col_map, c, r, m, t):
            noBrPointsPerEpoch = model.nbreakpoints
            nleaves = model.nleaves
            all_time_breakpoints, time_breakpoints = default_bps(model, c, r, t)

            M = []
            for e in xrange(len(noBrPointsPerEpoch)):
                newM = identity(nleaves)
                newM[:] = m[e]
                M.append(newM)

            pi, T, E = model.run(r, c, time_breakpoints, M, col_map=col_map)
            return pi, T, E

        pi, T, E = get_matrices(modelILS, COL_MAP, est_c, est_r, est_m, est_t)
        pi = array(pi, dtype=float64)
        T = array(T, dtype=float64)
        E = array(E, dtype=float64)
        E = E / E.sum(axis=1).reshape((E.shape[0],1))
        T = T / T.sum(axis=1)

        obs = all_obs[0]
        A, B, C1, logL = calc_forward_backward(pi,T,E, obs)

        winsize = int(self.winsize)
        L = int(self.L)

        def transitions_in_window(F, B, n, d, transitions, emissions, obs):
            m = n+d-1
            K = m-n
            Tnm = fill_table_fast(F, transitions, emissions, n, m, len(transitions), obs)
            #Tnm = fill_table(F, transitions, emissions, n, m, len(transitions), obs)
            res = zeros(K)
            for k in range(K):
                res[k] = (Tnm[:,k]*B).sum()
            #s = "\t".join(map(str, res))
            res = res / res.sum()
            #extra = ' '.join(map(str, around(res[0:10], 10)))
            res = res * arange(0, K)
            return res.sum()#, extra


        S = A.shape[1]
        B = B / B.sum(axis=1).reshape((L,1))

        for i, w in enumerate(range(0, L-winsize, winsize)):
            nrExpected = transitions_in_window(A[w,:], B[w+winsize-1,:], w, winsize, T, E, obs)
            print >>sys.stderr, nrExpected
            self.expectedTransitions.append(nrExpected)

# 
# import sys
# import os
# 
# import coalhmm
# from coalhmm.model import build_epoch_seperated_model
# from coalhmm.optimize import logL_multiseq, readObservations
# from scipy import *
# from scipy.optimize import fmin
# 
# # In this example the col map will be based on the actual observations in the file
# COL_MAP = dict()
# 
# 
# from coalhmm.mini_hmm import calc_forward
# def logLikelihood_for_Ronly(model, F, B, obs, c,r,m,t):
#     def ronly_forward(pi, T, E, obs):
#         forward_s, scales_s, logL_s = calc_forward(F, T, E, obs)
#         forward_s[-1,:] *= scales_s[-1]
#         forward_s[-1,:] = forward_s[-1,:] * B
#         scales_s[-1] = forward_s[-1,:].sum()
#         return log(scales_s).sum()
# 
#     def ronly_prepare(pi, T, E):
#         return (array(pi, dtype=float64),
#                 array(T, dtype=float64),
#                 array(E, dtype=float64))
#     return logL_multiseq(model, [obs], COL_MAP, c,r,m,t, prepare_matrices=ronly_prepare, single_logL=ronly_forward)
# 
# def optimize_f(model, obs, f_logL, init):
#     return fmin(lambda x: -f_logL(model, obs, *x), init, full_output=True, disp=False)
# 
# def estimate_Ronly(model, obs, F, B, T, C, R, outfile='/dev/stdout'):
#     def logL_all(model, obs, r):
#         r = exp(r)
#         if r <= 0:
#             return -1e18
#         res = logLikelihood_for_Ronly(model, F,B, obs, [C]*2, [r]*2, [0.0]*2, [0.0,T])
#         os.system("echo '%s' >> %s" % ('\t'.join(map(str, [r, res])), outfile))
#         return res
#     os.system("echo '%s' > %s" % ('\t'.join(map(str, ["R", "logL"])), outfile))
#     est, L, _, _, _ = optimize_f(model, obs, logL_all, (log(R),))
#     #       logL.   estimates
#     return (L, list([exp(est[0])]))
# 
# 
# for seq in seqs:
#     obs, _ = readObservations(seq, names, COL_MAP)
#     all_obs.append(obs)
# 
# L = int(1e6)
# winsize = int(100e3)
# stepsize = int(10e3)
# 
# obs = all_obs[0]
# 
# 
# S = len(modelI.tree_map)
# dummy = ones(S)/S
# 
# outfile = open("data/x_%i.logL" % (inst), 'w')
# maxL, est = estimate_Ronly(modelI, obs, dummy,dummy, i_t, i_c, i_r)
# estR = est[0]
# 
# for a in range(0, L-winsize, stepsize):
#     b = a + winsize
#     L1, est1 = estimate_Ronly(modelI, obs[a:b], dummy,dummy, i_t, i_c, i_r, outfile='/dev/null')
#     L2, est2 = estimate_Ronly(modelI, obs[a:b], dummy,dummy, i_t, i_c, estR, outfile='/dev/null')
#     print >>outfile, inst, a,b, L1,est1[0], L2,est2[0]
#     print            inst, a,b, L1,est1[0], L2,est2[0]
# 
# outfile.close()

sys.setrecursionlimit(100000)

class CoaSimSimulationHook(object):

    def __init__(self):
        self.recombinationTimes = []
        self.recombinationPoints = []
        self.recombinationLeaves = []
        self.trees = []

    def run(self, arg, trees, args):

        def collect_leaves(node):

            processed_nodes = set()

            def recurse(node):
                if isinstance(node, CoaSim.LeafNode):
                    return [node.sampleID]
                else:
                    res = []
                    for child in node.children:
                        if child.coreID not in processed_nodes:
                            res.extend(recurse(child))
                            processed_nodes.add(child.coreID)
                    return res

            return tuple(recurse(node))

        ######################################
        if coaSimLeafHack:
            # version of the below code commented out because this functionality has a memory leak in CoaSim
            # loop over nodes in the arg and collect recombination points and times:
            recombinationEvents = dict()
            for n in arg.nodes:            
                if isinstance(n, CoaSim.RecombinationNode) and n.isAncestral(n.recombinationPoint):
                    recombinationEvents[n.recombinationPoint] = (n.eventTime, sorted(collect_leaves(n)))

            # sort and scale points and times:
            for p, (t, l) in sorted(recombinationEvents.items()):
                self.recombinationPoints.append(int(round(p  * args["length"])))
                self.recombinationTimes.append(int(round(t * (2.0*args["NeRef"]) * args["g"])))
                self.recombinationLeaves.append(l)
        else:
            # loop over nodes in the arg and collect recombination points and times:
            recombinationEvents = dict()
            for n in arg.nodes:            
                if isinstance(n, CoaSim.RecombinationNode) and n.isAncestral(n.recombinationPoint):
                    recombinationEvents[n.recombinationPoint] = n.eventTime

            # sort and scale points and times:
            for p, t in sorted(recombinationEvents.items()):
                self.recombinationPoints.append(int(round(p  * args["length"])))
                self.recombinationTimes.append(int(round(t * (2.0*args["NeRef"]) * args["g"])))
        ######################################

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


class MaCSSimulationHook(object):

    def __init__(self):
        self.recombinationPoints = []
        self.trees = []

    def run(self, trees, args):
        for i, (s, e, t) in enumerate(trees):
            if i > 0:
                self.recombinationPoints.append(int(round(s  * args["length"])))
            self.trees.append(t)


# class ILSestimationHook(object):
# 
#     def __init__(self):
#         self.startCoord = []
#         self.endCoord = []
#         self.fromState = list()
#         self.toState = list()
#         self.infStates = list()
# 
#     def run(self, prefix, args):
#         prevInferredState = None
#         prevState = None
#         prevPos = None
#         f = open(prefix + ".posterior.tbl", 'r')
#         f.next() # remove header
#         prevPos = 0
#         pos = 0
#         for l in f:
#             lst = map(float, l.split()[1:5])
#             #argmax = lambda lst: max(izip(lst, xrange(len(lst))))[1]
#             maxProb = max(lst)
#             self.infStates.append(lst.index(maxProb))
#             if maxProb > 0.5:
#                 state = lst.index(maxProb)
#                 if prevState is not None and state != prevState:
#                     self.startCoord.append(prevPos)
#                     self.endCoord.append(pos)
#                     self.fromState.append(prevState)
#                     self.toState.append(state)
#                 prevState, prevPos = state, pos
#             pos += 1


class ILS09estimationHook(object):

    def __init__(self, minProb, maxSpan):
        self.minProb = minProb
        self.maxSpan = maxSpan
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
        f.next() # remove header
        prevPos = 0
        pos = 0
        ilsBases = 0
        non_ilsBases = 0

        for l in f:
            lst = map(float, l.split()[1:5])
            #argmax = lambda lst: max(izip(lst, xrange(len(lst))))[1]
            maxProb = max(lst)
            self.infStates.append(lst.index(maxProb))
            state = lst.index(maxProb)
            assert state <= 3, state
            if state <= 1:
                non_ilsBases += 1
            elif state >= 2:
                ilsBases += 1
                
            if maxProb > self.minProb:
                #if prevState is not None and state != prevState:
                if prevState is not None and state != prevState and pos - prevPos < self.maxSpan:
                    self.startCoord.append(prevPos)
                    self.endCoord.append(pos)
                    self.fromState.append(prevState)
                    self.toState.append(state)
                prevState, prevPos = state, pos
            pos += 1

        self.ilsBases = ilsBases
        self.non_ilsBases = non_ilsBases

# class CTMCestimationHook(object):
# 
#     def __init__(self):
#         self.startCoord = []
#         self.endCoord = []
#         self.fromState = list()
#         self.toState = list()
#         self.infStates = list()
# 
#     def run(self, prefix, args):
#         f = open("%s_model.pickle")
#         model = pickle.load(f)
#         f.close()
#         tree_map = model.treemap
#         # PARSE MODEL TO SO WE CAN DISTINGUISH MODELS WITH DIFFERENT TOPOLOGIES AND WIGHIN THOSE WHAT THE BRANCHLENGTHS (IN TERMS OF BINS) ARE....
# 
#         prevInferredState = None
#         prevState = None
#         prevPos = None
#         f = open(prefix + "_pds.txt", 'r')
#         f.next() # remove header
#         prevPos = 0
#         pos = 0
#         for l in f:
#             lst = map(float, l.split()[1:5])
#             #argmax = lambda lst: max(izip(lst, xrange(len(lst))))[1]
#             maxProb = max(lst)
#             self.infStates.append(lst.index(maxProb))
#             if maxProb > 2.0 / args["nstates"]:
#                 state = lst.index(maxProb)
#                 if prevState is not None and state != prevState:
#                     self.startCoord.append(prevPos)
#                     self.endCoord.append(pos)
#                     self.fromState.append(prevState)
#                     self.toState.append(state)
#                 prevState, prevPos = state, pos
#             pos += 1

class TreeHeight(tree.TreeVisitor):
    def __init__(self):
        self.height = 0
        self.max = 0
        self.depths = list()
        
    def reset(self):
        TreeHeight.__init__(self)
# 
#     def result(self):
#         return self.max
# 
#     def get_result(self):
#         return self.max

    def pre_visit_edge(self,src,bootstrap,length,dst):
        self.height = self.height + length
        self.max = max(self.max, self.height)
        
    def post_visit_edge(self,src,bootstrap,length,dst):
        self.height = self.height - length
        self.depths.append(self.height)
        
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
        th.depths = []
        lc.count = 0
        t.dfs_traverse(th)
        t.dfs_traverse(lc)

#         print t
#         print [th.max - x for x in th.depths]

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
    assert abs(1 - mapList[-1][0]) < 1e-10, 1 - mapList[-1][0]
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

    for i, (interval_start, interval_end, tree) in enumerate(intervals):
#     for i, interval in enumerate(intervals):
#         interval_start, interval_end = interval.start, interval.end
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
#         newIntervals.append((start, end, parse_tree(str(interval.tree))))
        newIntervals.append((start, end, parse_tree(str(tree))))

    # this handles the special case where the last specified recombination rate segment
    # has rate zero:
    if mapList[-1][2] == 0:
        newIntervals[-1] = (newIntervals[-1][0], newIntervals[-1][1] + mapList[-1][1], newIntervals[-1][2])

#     with open("test.tbl", 'w') as f:
#         for start, end, tree in newIntervals:
#             print >>f, start + (end-start)/2.0
#     sys.exit()
    
    return newIntervals


def reRootMaCSTree(tree, **args):
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
    val = args["tau123"] - max
    s = '(' + s + ' : ' + str(val) + " , '3' : " + str(args["tau123"]) + ' );'
    return s

def simulateMaCS(**args):

#     args["tau1"] = args["T1"] / (2*args["NeRef"])
#     args["tau12"] = args["T12"] / (2*args["NeRef"])
#     args["tau123"] = args["T123"] / (2*args["NeRef"])

    args["tau1"] = args["T1"] / args["g"] / (2*args["NeRef"])
    args["tau12"] = args["T12"] / args["g"] / (2*args["NeRef"])
    args["tau123"] = args["T123"] / args["g"] / (2*args["NeRef"])


    outputdir = tempfile.mkdtemp()

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)



#         rescale <- function(data, subst.per.year, years.per.gen) {
#   subst.per.gen <- years.per.gen*subst.per.year
# 
#   data$T1 <- data$tau1 / subst.per.year
#   data$T12 <- (data$tau1+data$tau2) / subst.per.year
#   data$T123 <- (data$tau1+data$tau2+data$tau3) / subst.per.year
# 
#   data$Ne12 <- data$theta1 / (2*subst.per.gen)
#   data$Ne123 <- data$theta2 / (2*subst.per.gen)
#   data$rec.rate <- data$rho * subst.per.gen * 1e6 * 1e2 # cM/Mb
#   return(data)
# }
# 
# "tau1=" + str(T1 * u * 2 * NeRef * g) + " tau2=" + str((T12-T1) * u * 2 * NeRef * g) + " c2=" +\
#         str(((T123- T12) * NeRef - 4/3 * N123) * u * 2 * g) + " theta1=" + str(N12 * u * 2 * g) + " theta2=" + \
#         str(N123 * u * 2 * g) + " rho=" + str(r / (u * g)) + " theta12=" + str(N12 * u * 2 * g) + " DATA=" + sequence
# 


    cmd = macs_exe + " 3 " + str(args["length"]) + " -T -r " + str(4 * args["NeRef"] * args["r"]) + " -I 3 1 1 1 " + \
    " -n 1 " + str(args["N1"] / args["NeRef"]) + " -n 2 " + str(args["N1"] / args["NeRef"]) + " -n 3 " + str(args["N1"] / args["NeRef"]) +  \
    " -ej " + str(args["tau1"]) + " 2 3 -en " + str(args["tau1"]) + " 2 " + str(args["N12"] / args["NeRef"]) +  " -en " + str(args["tau1"]) + " 3 " + str(args["N12"] / args["NeRef"]) + \
    " -ej " + str(args["tau12"] ) +" 3 1 -en " + str(args["tau12"]) + " 1 " + str(args["N123"] / args["NeRef"]) 
#    print cmd

    p = subprocess.Popen(cmd, env=os.environ, shell=True, cwd=macs_dir, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    #p = subprocess.Popen(cmd, env=os.environ, shell=True, cwd=macs_dir, stdout=subprocess.PIPE)
    stdout, stderr = p.communicate()
    stderr = "" # del stderr

    trees = list()

    treelengthmatcher = re.compile(': *-?\d*\.?\d*e?-?\d*')
    macsmatch = re.compile('.*\[(.*?)\](.*)')
    #macsmatch = re.compile('.*\[(\d*)\](.*)')

#    TEST = list()
    prevTree = None
    linecount = 0
#    tmptree = open("tmptree" , "w")
    last = 0
#    for item in fin :
    for item in stdout.strip().split('\n') :

        linecount += 1
        #first 2 lines does not contain trees
        if(linecount < 3):
            continue

        m = macsmatch.match(item)
#        l = int(m.group(1))
        l = float(m.group(1))
        tree = m.group(2)
        tree = tree.replace("0:", "'0':")
        tree = tree.replace("1:", "'1':")
        tree = tree.replace("2:", "'2':")
        
        # add outgroup
        if "addOutGroup" in args and args["addOutGroup"]:
            tree = reRootMaCSTree(tree, **args)
        
        start = float(last)/args["length"]
        end = float(last+l)/args["length"]



        # scale trees
        lens = treelengthmatcher.findall(tree)
        for le in lens :
            fl = float(le.replace(':', '')) * 2 * args["NeRef"] * args["u"] * args["g"]
            tree = tree.replace(le, ':' + str(fl))
        
        last = last + l

        # since we are not anaylzing the recombinations we merge intervals with
        # identical trees - correponds to removing dimond recombinations one external
        # branches.
        if tree == prevTree:
            trees[-1][1] = end
        else:
            trees.append([start, end, tree])
        prevTree = tree

    # HACK: even with my patch MaCS does not produce intervals up to exactly one -
    # though very nearly so. So we help it along for the last bit:
    trees[-1][1] = 1.0

#     # HACK: MaCS makes a rounding error with rounding to an int base position. All the
#     # ones that should be rounded up are rounded down. So the total nr of bases is about
#     # nr.rec.events/2 short. We fix that here:
#     factor = 1/float(trees[-1][1])
#     for i in range(len(trees)):
#         trees[i][0] = trees[i][0]*factor
#         trees[i][1] = trees[i][1]*factor
        

    if "recMap" in args:
        # map intervals and parse trees
        trees = mapIntervalsAndParseTrees(trees, args["recMap"])
    else:
        trees = [(s, e, parse_tree(t)) for s, e, t in trees]

    if "hook" in args:
        args["hook"].run(trees, args)

#     ####
#     def getTree(tree):
#         left, right = tree.get_edges()
#         if isinstance(left[0], Leaf) and left[0].identifier == '3': # mayby identifier is not an int but a str...
#             # outgroup is included - call on ingroup
#             return right[0]
#         elif isinstance(right[0], Leaf) and right[0].identifier == '3': # mayby identifier is not an int but a str...
#             # outgroup is included - call on ingroup
#             return left[0]
#         else:
#             return tree
# 
#     th = TreeHeight()
# 
#     for (start, end, tree) in trees:
#         th.reset()        
#         t = getTree(tree)
#         t.dfs_traverse(th)
#         h = th.get_result()
#         print h
# 
#     with open("macs.trees", 'w') as f:
#         for s, e, t in trees:
#             print >>f, s, e, t
#    sys.exit()
#    ####

    return bppSeqGen(trees, **args)


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

     
    #run the stuff
    
    st = spec(**args)
    kei = True
    if "keepEmptyIntervals" in args:
        kei = args["keepEmptyIntervals"];
    
#    print "simulating"
        
    arg = CoaSim.simulate([], st, rho = 4 * args["NeRef"] * args["length"] * args["r"],
            keepEmptyIntervals=kei, migrationSpec=args["migration"]) #what about 4

    ## # FIXME: possibly do this instead: collapse all runs of identical trees:
    ## trees = list()
    ## prevTree = None            
    ## for tree in [(x.start, x.end, str(x.tree)) for x in arg.intervals]:
    ##     if tree == prevTree:
    ##         trees[-1][1] = end
    ##     else:
    ##         trees.append([start, end, tree])
    ##     prevTree = tree

    ## if "recMap" in args:
    ##     # map intervals and parse trees
    ##     trees = mapIntervalsAndParseTrees(trees, args["recMap"])
    ## else:
    ##     # parse trees
    ##     for i in range(len(trees)):
    ##         trees[i][2] = parse_tree(trees[i][2])        

    
    if "recMap" in args:
        # map intervals and parse trees
        trees = mapIntervalsAndParseTrees([(x.start, x.end, x.tree) for x in arg.intervals], args["recMap"])
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

#     for t in trees:
#         print t

# ####
#     def getTree(tree):
#         left, right = tree.get_edges()
#         if isinstance(left[0], Leaf) and left[0].identifier == '3': # mayby identifier is not an int but a str...
#             # outgroup is included - call on ingroup
#             return right[0]
#         elif isinstance(right[0], Leaf) and right[0].identifier == '3': # mayby identifier is not an int but a str...
#             # outgroup is included - call on ingroup
#             return left[0]
#         else:
#             return tree
# 
#     th = TreeHeight()
# 
#     for (start, end, tree) in trees:
#         th.reset()        
#         t = getTree(tree)
#         t.dfs_traverse(th)
#         h = th.get_result()
#         print h
# 
# 
#     with open("coasim.trees", 'w') as f:
#         for s, e, t in trees:
#             print >>f, s, e, t
#    sys.exit()
####


    return bppSeqGen(trees, **args)





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
#     sys.exit()
    p = subprocess.Popen(cmd + " " + param, env=os.environ, shell=True, cwd=coalhmm_dir, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    stdout, stderr = p.communicate()
    sleep(30)
    #print stdout
    #print stderr
    if "Exception" in stderr or "Oups... something abnormal happened!" in stderr or "Optimization failed because likelihood function returned NaN" in stderr:
        print stderr
        return None
    if "Exception" in stdout:
        print stdout
        return None
#    sequence = '/var/folders/8j/_27xl8wd6vn5krws659zcr1h0000gn/T/tmpAkZskf'
    if "hook" in args:
        args["hook"].run(sequence, args)
    
    input_t = Table.Table().add_row(_tau1=T1 * u * 2 * NeRef * g, _tau2=(T12-T1) * u * 2 * NeRef * g, \
                              _c2=((T123- T12) * NeRef - 4/3 * N123) * u * 2 * g, _theta1=N12 * u * 2 * g, \
                              _theta2=N123 * u * 2 * g, _rho=r / (u * g), _theta12=N12 * u * 2 * g)

    return Table.Table().load_kv(sequence + ".user.txt").join(input_t)


def estimate_ctmcILS(sequence, **args):

    likelihood, estimates, startvalues = runILSctmc([sequence], **args)

    u, g = args['u'], args['g']
    return Table.Table().add_row(_u = u,
                                 _g = g,
                                 _r = startvalues["i_r"]*u,
                                 _c1 = startvalues["i_c1"],
                                 _c2 = startvalues["i_c2"],
                                 _c3 = startvalues["i_c3"],
                                 _N1 = 1.0/(2*g*u*startvalues['i_c1']),
                                 _N12 = 1.0/(2*g*u*startvalues['i_c2']),
                                 _N123 = 1.0/(2*g*u*startvalues['i_c3']),
                                 _T12 = startvalues['i_t1'] * 2 * g / u,
                                 _T123 = startvalues['i_t2'] * 2 * g / u,
                                 r = estimates["i_r"]*u,
                                 c1 = estimates["i_c1"],
                                 c2 = estimates["i_c2"],
                                 c3 = estimates["i_c3"],
                                 N1 = 1.0/(2*g*u*estimates['i_c1']),
                                 N12 = 1.0/(2*g*u*estimates['i_c2']),
                                 N123 = 1.0/(2*g*u*estimates['i_c3']),
                                 T12 = estimates['i_t1'] * 2 * g / u,
                                 T123 = estimates['i_t2'] * 2 * g / u)

