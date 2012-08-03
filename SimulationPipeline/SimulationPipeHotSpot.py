import CoaSim, os, tempfile, table, sys
from CoaSim.popStructure import Population as P, Sample as S, Merge as M
from newick import *

coalhmm_exe = "~tth/local/bin/coalhmm --noninteractive=yes"
bppseqgen_exe = "~tth/local/bin/bppseqgen --noninteractive=yes"

# coalhmm_exe = "/usr/local/bin/coalhmm --noninteractive=yes"
# bppseqgen_exe = "/usr/local/bin/bppseqgen --noninteractive=yes"

# coalhmm_exe = "coalhmm --noninteractive=yes"
# bppseqgen_exe = "bppseqgen --noninteractive=yes"

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

    print path1

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

    print cmd
    os.system(cmd)
    #return path to sequence
    #todo: remove tmp file
    return path2

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
    

    trees = []
    if "hotspotWidth" in args and "hotspotRatio" in args:
        # hotspotWidth
        hw = args["hotspotWidth"]
        # ratio of hotspot non-hotstop recombination rate:
        F = args["hotspotRatio"]
        p = hw / (1 - hw)
        # width of hotspot region before contraction:
        ohw = (F*p)/(F*p+1)
        # expansion factor:
        e = (1-hw) / (1-ohw)
        # contraction factor:
        c = hw / ohw

        args["hotSpotBackgroundRate"] = args["r"] * (1-hw)/(1-ohw)
        print >>sys.stderr, "Simulating hotspot", args["hotSpotBackgroundRate"]

#         # nr hotspots:
#         nrHotSpots = 4
#         # o iohw into 
#         m = generateMap(nrHotSpots, ohw, c, e):
        
        # map
        m = [(0, e), ((1-ohw)/2.0, c), (((1-ohw)/2.0)+ohw, e), (1.0, None)]
        j = 0    
        prevEnd = 0
        recPoints = list()
        for i, interval in enumerate(arg.intervals):
            while j < len(m)-1 and interval.start >= m[j+1][0]:
                j += 1
            start = prevEnd
            if j < len(m)-1 and interval.end > m[j+1][0]:
                end = prevEnd + (m[j+1][0] - interval.start) * m[j][1] + (interval.end - m[j+1][0]) * m[j+1][1]
            else:
                end = prevEnd + (interval.end - interval.start) * m[j][1]
            prevEnd = end
            trees.append((start, end, parse_tree(str(interval.tree))))
        
    else:
        for interval in arg.intervals:
            trees.append((interval.start, interval.end, parse_tree(str(interval.tree))))        
#     for s, e, t in trees:
#         print s, e
#     import sys
#     sys.exit()

    #reroot trees
    if "addOutGroup" in args and args["addOutGroup"]:
        trees = addOutGroup(trees, **args)
    

#    print "scaling"
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
        args["hook"].parseARG(arg, trees, args)
#         args["hook"].parseTrees(arg, args)
    
    return bppSeqGen(trees, **args)

def estimate(sequence, **args):
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

    print cmd + " " + param
    os.system(cmd + " " + param)

    if "hook" in args:
        args["hook"].run(sequence, args)

    input_t = table.Table().add_row(_tau1=T1 * u * 2 * NeRef * g, _tau2=(T12-T1) * u * 2 * NeRef * g, \
                              _c2=((T123- T12) * NeRef - 4/3 * N123) * u * 2 * g, _theta1=N12 * u * 2 * g, \
                              _theta2=N123 * u * 2 * g, _rho=r / (u * g), _theta12=N12 * u * 2 * g)

    return table.Table().load_kv(sequence + ".user.txt").join(input_t)
    
'''
    T = []
    t = "T"
    for i in range (1, count):
        t = t + str(i)
        tmp1 = args[t]
        tmp2 = T.pop();
        T.append(tmp2)
        T.append(tmp1 - tmp2)
'''

# #begin new model
# 
# from scipy import *
# from scipy.linalg import expm
# 
# #from time_plot import *
# 
# import sys
# sys.path.append( "nymodel" )
# sys.path.append( "/users/asand/public_html/" )
# 
# from model import build_epoch_seperated_model, build_simple_model
# from optimize import logLikelihood, readObservations
# from scipy.optimize import fmin
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
# def optimize_f(file_in, model, seqnames, f_logL, init):
#     obs, colmap = readObservations(file_in, seqnames)
#     nbps = model.nbreakpoints
#     return fmin(lambda x: -f_logL(model, obs, colmap, *x), init)
# 
# def estimate_untitled_3_all(seq, T1, T2, C, R, nstates=5):
#     def logL_all(model, obs, colmap, c):
#         if min([c]) <= 0:
#             return -1e18
# 	print "C =", c
#         return logLikelihood(model, obs, colmap, [c,c], [R,R], [0.0,0.0], [0.0,T1])
# 
#     #antal epoker i model 2, hvad de konvergerer til 0, antal states i de epoker hhv. 1, nstates
#     ssplit = build_epoch_seperated_model(2, [[0,0]], [1,nstates], None)
#     names = ["'0'", "'1'"]
#     est = optimize_f(seq, ssplit, names, logL_all, (C,))
#     
#     obs, colmap = readObservations(seq, names)
#     _, pd, tbps, pi = logLikelihood(ssplit, obs, colmap, [est[0],est[0]], [R,R], [0.0,0.0], [0.0,T1], posterior_decoding=True)
# 
#     print pd.get_no_rows(), "x", pd.get_no_columns()
# 
#     f_bps = open("bps.txt", "w")
#     print >>f_bps, "\t".join(map(str, tbps))
#     f_bps.close()
#     f_pi = open("pi.txt", "w")
#     print >>f_pi, "\t".join([str(pi[x]) for x in xrange(len(pi))])
#     f_pi.close()
# 
#     pd_file = open("pds.txt", "w")
#     print >>pd_file, matrix_to_string(pd)
#     pd_file.close()
# 
#     return table.Table().add_row(C=est[0],R=R,tau_1=T1)
