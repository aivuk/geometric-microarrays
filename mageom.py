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
    n_genes = len(g_exp)
    e_g = np.zeros(n_genes)
    comb_c_1 = list(i.combinations(classes[0], 2))
    comb_c_2 = list(i.combinations(classes[1], 2))
    prod_f = list(i.product(classes[0], classes[1]))
    len_c = len(comb_c_1) + len(comb_c_2)
    len_f = len(prod_f)
    
    e_g = np.zeros(n_genes)
    for g in range(0, n_genes):
        e_g[g] = 1./len_c*(np.sum(map(lambda (i,j): np.sqrt((g_exp[g][i] - g_exp[g][j])**2), comb_c_1)) + \
                            np.sum(map(lambda (i,j): np.sqrt((g_exp[g][i] - g_exp[g][j])**2), comb_c_2))) \
               - 1./len_f*np.sum(map(lambda (i,j): np.sqrt((g_exp[g][i] - g_exp[g][j])**2), prod_f)) 
    
    return np.argsort(e_g)
