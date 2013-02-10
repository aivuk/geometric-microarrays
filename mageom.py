import numpy as np
import itertools as i

def calc_metric_tensor(g_exp, classes):
    """
    @type g_exp: double dimensional matrix
    @param g_exp: one experiment per line and one gene per column
    @type classes: list with exactly two lists
    @param classes: a list with two element, first a list with 
    @rtype: one dimensional vector 
    @return: the metric tensor
    """
    comb_c = np.array(list(i.combinations(classes[0], 2)) + list(i.combinations(classes[1], 2)))
    prod_f = np.array(list(i.product(classes[0], classes[1])))
    
    e_g = 1./len(comb_c)*np.sum(np.abs(g_exp[comb_c[:,0],:] - g_exp[comb_c[:,1],:]), axis=0) \
        - 1./len(prod_f)*np.sum(np.abs(g_exp[prod_f[:,0],:] - g_exp[prod_f[:,1],:]), axis=0) 

    return e_g
