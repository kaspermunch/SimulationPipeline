

import sys
import os

import coalhmm
from coalhmm.model import build_epoch_seperated_model
from coalhmm.optimize import logL_multiseq, readObservations
from scipy import *
from scipy.optimize import fmin

from pyZipHMM import *
from itertools import product
import tempfile
from coalhmm.optimize import default_bps
from coalhmm.mini_hmm import calc_forward_backward

def transitions(T, B, S, K):
    assert T.shape == (S,K)
    assert B.shape == (S,)
    assert S > 0 and K > 0

    result = zeros(K)
    for k in range(K):
        result[k] = (T[:,k]*B).sum()
    result = result / result.sum()
    result = result * arange(1, K+1)
    return result.sum()

def read_matrix(fname, mapping=float):
    f = open(fname)
    res = []
    for line in f:
        vals = list(map(mapping, line.split()))
        res.append(vals)
    f.close()
    res = array(res)
#    print "%20s" % fname, res.shape
    return res


COL_MAP = dict()


def zipHMM_prepare_matrices(mpi, mT, mE):
#    print '.',
    sys.stdout.flush()
    pi = Matrix(mpi.shape[0], 1)
    for i in range(mpi.shape[0]):
        pi[0, i] = mpi[i]
    T  = Matrix(mT.shape[0], mT.shape[1])
    for i in range(mT.shape[0]):
        for j in range(mT.shape[1]):
            T[i, j] = mT[i, j]
    E  = Matrix(mE.shape[0], mE.shape[1])
    for i in range(mE.shape[0]):
        for j in range(mE.shape[1]):
            E[i, j] = mE[i, j]
    return pi, T, E

def zipHMM_single_logL(pi, T, E, obs):
    return obs.forward(pi, T, E)

def logLikelihood(model, all_obs, c,r,m,t):
    return logL_multiseq(model, all_obs, COL_MAP, c,r,m,t,
            prepare_matrices=zipHMM_prepare_matrices,
            single_logL=zipHMM_single_logL)

def optimize_f(model, obs, f_logL, init):
    return fmin(lambda x: -f_logL(model, obs, *x), init, full_output=True)

def estimate_ILS(model, obs, T1, T2, C1, C2, C3, R, outfile="/dev/null"):

    def logL_all(model, all_obs, t1, t2, c1, c2, c3, r):
        if min([t1,t2,c1,c2,c3,r]) <= 0 or t2 <= t1:
            return -1e18
        res = logLikelihood(model, obs, [c1, c2, c3], [r]*3, [0.0]*3, [0.0,t1,t2])
        os.system("echo '%s' >> %s" % ('\t'.join(map(str, [t1,t2,c1,c2,c3, r, res])), outfile))
        return res

    os.system("echo '%s' > %s" % ('\t'.join(map(str, ["t1","t2", "c1", "c2", "c3", "r", "logL"])), outfile))
    est, L, _, _, _ = optimize_f(model, obs, logL_all, (T1,T2,C1,C2,C3,R))
    #       logL.   estimates
    return (L, list(est))

def estimate_I(model, obs, T, C, R, outfile="/dev/null"):
    def logL_all(model, all_obs, t, c, r):
        if min([t,c,r]) <= 0:
            return -1e18
        res = logLikelihood(model, obs, [c]*2, [r]*2, [0.0]*2, [0.0,t])
        os.system("echo '%s' >> %s" % ('\t'.join(map(str, [t,c, r, res])), outfile))
        return res
    os.system("echo '%s' > %s" % ('\t'.join(map(str, ["t", "c", "r", "logL"])), outfile))
    est, L, _, _, _ = optimize_f(model, obs, logL_all, (T,C,R))
    #       logL.   estimates
    return (L, list(est))

def runILSctmc(seqs, **args):
    global COL_MAP
    NeRef = args["NeRef"]
    g = args["g"]
    u = args["u"]
    i_r = args["r"]
    i_N1 = args["N1"] 
    i_N2 = args["N12"]
    i_N3 = args["N123"]
    i_c1 = 1.0/(2*g*u*i_N1)
    i_c2 = 1.0/(2*g*u*i_N2)
    i_c3 = 1.0/(2*g*u*i_N3)
    i_t1 = u * args["T12"] / g / (2*NeRef)
    i_t2 = u * args["T123"] / g / (2*NeRef)

    startvalues = dict( [ (name,eval(name)) for name in ['i_t1', 'i_t2', 'i_c1', 'i_c2', 'i_c3', 'i_r'] ] )

#     # generation time and mutation rate
#     g = 20
#     u = 1e-9
# 
#     # initial values for recombination rate, population size/coalescence rate and
#     # migration rate.
#     i_r = 0.4
#     i_N = 100e3
#     i_c = 1.0/(2*g*u*i_N)
#     i_t1 = 3.7e6*u
#     i_t2 = 5.95e6*u

    estates = 4
    print('running with %d epoch states' % estates, file=sys.stderr)
    
#     inst = int(sys.argv[1])
# 
#     seqs = ["data/x_ils_%i.fa" % inst]

    modelILS = build_epoch_seperated_model(3, [[0,0,1], [0,0]], [1,estates,estates])
    nstates = len(modelILS.tree_map)
    names = ["'1'", "'2'", "'0'"]
    #names = ["'0'", "'2'", "'1'"]

    all_obs = []
    forwarders = []

    COL_MAP = dict((v,i) for i,v in enumerate(product('ACGT', repeat=3)))
    for seq in seqs:
        # each sequence creates a new column map extended with any new columns that
        # might be seen.
        obs, colmap = readObservations(seq, names, COL_MAP)
        all_obs.append(obs)

#     print len(COL_MAP)

    doEstimate = False
    
    L = None
    estimates = list()
    if doEstimate:
        for obs in all_obs:
            print('next obs:')
            ffd, foutname = tempfile.mkstemp()
            print('  temp fd/name:', ffd, foutname)
            fout = os.fdopen(ffd, 'w')
            L = len(obs)
            for j in range(L-1):
                o = obs[j]
                print(o, end=' ', file=fout)
            print(obs[j], file=fout)
            fout.close()
            print('  written, creating forwarder')
            f = Forwarder.fromSequence(seqFilename = foutname, alphabetSize = len(COL_MAP))
            #f = Forwarder(seqFilename = foutname, nStates = len(modelILS.tree_map), nObservables = len(COL_MAP))
            print('  - done.')
            forwarders.append(f)
            os.system("rm %s" % foutname)                                                                

        L, estimates = estimate_ILS(modelILS, forwarders, i_t1, i_t2, i_c1, i_c2, i_c3, i_r, outfile="/dev/null")
        i_t1, i_t2, i_c1, i_c2, i_c3, i_r = estimates

    estimates = dict( [ (name,eval(name)) for name in ['i_t1', 'i_t2', 'i_c1', 'i_c2', 'i_c3', 'i_r'] ] )

    #print 'Estimates:'
    #print "\t".join(map(str, [L] + est))


    if "hook" in args:
        args["hook"].run(modelILS, COL_MAP, all_obs, [i_c1, i_c2, i_c3], [i_r]*3, [0]*3, [0, i_t1, i_t2])

    return L, estimates, startvalues
