import os
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import hicstuff.view as hcv
import hicstuff.io as hio


"""
!!! Warning build hicstuff from source to 
     use with Python 3.10 and newer !!!
"""

if __name__ == "__main__":

    data_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/data/'
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']
    outputs_dir = data_dir + 'outputs/'
    sparse_dir = outputs_dir + "sparse/"
    hic_dense_dir = outputs_dir + "dense/"

    to_dense = lambda x: hcv.sparse_to_dense(x, remove_diag=False)

    if not os.path.exists(hic_dense_dir):
        os.makedirs(hic_dense_dir)

    for sshic_dir in sshic_pcrdupt_dir:
        samples_dir = sparse_dir + sshic_dir
        samples = np.unique(os.listdir(samples_dir))

        for samp in samples:
            sparse_mat = hio.load_sparse_matrix(samples_dir+samp)
            dense_mat = to_dense(sparse_mat)