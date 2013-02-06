import numpy as np
import itertools as i
cimport numpy as np
cimport cython


@cython.boundscheck(False)
def calc_metric_tensor(g_exp, classes):
    """
    @type g_exp: double dimensional matrix
    @param g_exp: one experiment per line and one gene per column
    @type classes: list with exactly two lists
    @param classes: a list with two element, first a list with 
    @rtype: one dimensional vector 
    @return: the metric tensor
    """
    cdef int n_genes = g_exp.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] e_g = np.zeros(n_genes)
    comb_c_1 = list(i.combinations(classes[0], 2))
    comb_c_2 = list(i.combinations(classes[1], 2))
    prod_f = list(i.product(classes[0], classes[1]))
    cdef int len_c = len(comb_c_1) + len(comb_c_2)
    cdef int len_f = len(prod_f)
    cdef int g 
    for g in range(0, n_genes):
        e_g[g] = 1./len_c*(np.sum(map(lambda d: np.sqrt((g_exp[d[0]][g] - g_exp[d[1]][g])**2), comb_c_1)) + \
                            np.sum(map(lambda d: np.sqrt((g_exp[d[0]][g] - g_exp[d[1]][g])**2), comb_c_2))) \
               - 1./len_f*np.sum(map(lambda d: np.sqrt((g_exp[d[0]][g] - g_exp[d[1]][g])**2), prod_f)) 
    
    return e_g

def test_speed():
    cdef np.ndarray[np.float64_t, ndim=2] x = np.random.randn(200, 23000)
    c1 = range(0,100)
    c2 = range(100,200)
    cs = [c1, c2]
    return calc_metric_tensor(x, cs)

