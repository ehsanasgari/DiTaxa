#!/usr/bin/env python

from svm import *
from svmutil import *
import croc
import argparse
import utils

import random
random.seed(1982)


def ROC( ty, pv ):
    if len(ty) != len(pv):
        raise ValueError("len(ty) must equal to len(pv)")
    SD = croc.ScoredData()
    for t,p in zip(ty,pv):
        #print t,p
        SD.add( p,t )
    it = SD.sweep_threshold()
    #print it
    curve = croc.ROC( it )
    return (curve.area(),curve)


def get_best_param( res_cv ):
    cmax, gmax, rmax = None, None, -1.0
    for c,gdict in res_cv.items():
        for g,r in gdict.items():
            if rmax < r:
                rmax = r
                cmax = c
                gmax = g
    return cmax,gmax,rmax

def ms_cv( prob, param, nfolds = 10, metric = 'auc', c_list = None, g_list = None ):

    prob = zip(*prob)
    folds = get_folds( prob, nfolds )
    labels,vals,targets = [],[],[]
    c_orig, g_orig = param.C, param.gamma

    for nf,fold in enumerate(folds):
        training_folds = []
        for j,f in enumerate(folds):
            if nf != j:
                training_folds += f
    
        y, x = zip(*training_folds)
   
        c_list = c_list if c_list else [c_orig]
        g_list = g_list if g_list else [g_orig]
        
        cv_res = dict([(c,{}) for c in c_list])
        for c in c_list:
            for g in g_list:
                param.C = c
                param.gamma = g
                cv_res[c][g] = cv( (y, x), param, nfolds = nfolds, metric = metric )[0] 

        param.C, param.gamma, rmax = get_best_param(cv_res) 
        #sys.stderr.write( str(param.C) + str(param.gamma) + str(rmax) )

        m = svm_train( svm_problem( y, x ), param )
        yp, xp = zip(*fold)
        labs,acc,vs = svm_predict( yp, xp, m, options="-b 1" if param.probability else "" ) 
        
        labels += labs
        vals += vs
        targets += yp
    
    if metric == 'acc':
        ACC, MSE, SCC = evaluations(labels, targets)
        return ACC, None
    if metric == 'auc':
        auc,curve = ROC( targets, [-v[0] for v in vals]) 
        return auc,curve


def get_folds( prob, nfolds ):
    
    folds = [[] for i in range(nfolds)]
    random.shuffle( list(prob) )
    prob = sorted( prob, key = lambda x:x[0])
    
    for i,v in enumerate( prob ):
        folds[i%nfolds].append(v)
    
    return folds


def cv( prob, param, nfolds = 10, metric = 'acc' ):
    assert len(prob[0]) > nfolds

    prob = zip(*prob)

    folds = get_folds( prob, nfolds )
  
    labels,vals,targets = [],[],[]

    for nf,fold in enumerate(folds):
        training_folds = []
        for j,f in enumerate(folds):
            if nf != j:
                training_folds += f
       
        y, x = zip(*training_folds)
        m = svm_train( svm_problem( y, x ), param )
        yp, xp = zip(*fold)
        labs,acc,vs = svm_predict( yp, xp, m, options="-b 1" if param.probability else "" ) 
        labels += labs
        vals += vs
        targets += yp
    
    if metric == 'acc':
        ACC, MSE, SCC = evaluations(labels, targets)
        print "Acc",ACC
        return ACC,None
    if metric == 'auc':
        auc,curve = ROC( targets, [-v[0] for v in vals]) 
        print "AUC",auc
        return auc,curve
    

def read_params(args):
    parser = argparse.ArgumentParser(description='LibSVM extensions')
    arg = parser.add_argument
    arg( 'inp_f', metavar='INPUT_FILE', default=None, type=str )
    arg( '--libsvm', default=[], nargs='+', type =str )
    arg( '-n', default=10, type = int )
    arg( '-m', default='acc', type = str )
    arg( '-c', default=None, nargs='+', type =str )
    arg( '-g', default=None, nargs='+', type =str )
    arg( '--roc', default=None, type =str )
    return vars(parser.parse_args())
        
args = read_params(sys.argv)

p = svm_read_problem(args['inp_f'])
print " ".join(args['libsvm']).replace("_","-")
param = svm_parameter( " ".join(args['libsvm']).replace("_","-")  )


c_list = [float(f) for f in args['c']] if args['c'] else [1.0]
g_list = [float(f) for f in args['g']] if args['g'] else [1.0]

v,g = ms_cv( p, param, nfolds = args['n'], metric = args['m'], c_list = c_list, g_list = g_list )

if args['m'] == 'auc':
    print "CV AUC\t"+str(v)
elif args['m'] == 'acc':
    print "CV acc\t"+str(v)

if args['roc'] and args['m'] == 'auc':
    with open( args['roc'], "w" ) as out:
        for v in g:
            out.write( str(v[0]) + "\t" + str(v[1]) + "\n" )


