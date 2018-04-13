__author__ = "Ehsaneddin Asgari"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu or ehsaneddin.asgari@helmholtz-hzi.de"
__project__ = "LLP - Life Language Processing"
__website__ = "https://llp.berkeley.edu/"

from scipy import stats
from sklearn.preprocessing import normalize
import numpy as np
from sklearn.preprocessing import  normalize


def get_sym_kl_rows(A):
    '''
    :param A: matrix A
    :return: Efficient implementation to calculate kl-divergence between rows in A
    '''
    norm_A=normalize(A+np.finfo(np.float64).eps, norm='l1')
    a=stats.entropy(norm_A.T[:,:,None], norm_A.T[:,None,:])
    return a+a.T

def get_kl_rows(A):
    '''
    :param A: matrix A
    :return: Efficient implementation to calculate kl-divergence between rows in A
    '''
    norm_A=normalize(A+np.finfo(np.float64).eps, norm='l1')
    return stats.entropy(norm_A.T[:,:,None], norm_A.T[:,None,:])

def normalize_mat(A,norm ='l1', axis=1):
    '''

    :param A:
    :param norm:
    :param axis: 0 colum
    :return:
    '''

    return normalize(A, norm=norm, axis=axis)
