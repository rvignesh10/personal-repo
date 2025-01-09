import numpy as np
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt

def get_csr_case1():
    # return car info: case 1
    nrows = 3
    ncols = nrows

    row = np.array([0, 1, 2, 0, 1, 0])
    col = np.array([0, 1, 2, 1, 1, 1])
  
    assert(np.max(row)<=nrows-1)
    assert(np.max(col)<=ncols-1)
    assert(len(row)==len(col))

    # data with ones (not zeros)
    data = np.ones(len(row))

    return nrows, ncols, row, col, data

def get_row_col(ien):
    # return row and col for given ien

    Ne = ien.shape[0]
    nen = ien.shape[1]

    row = np.zeros(int(Ne*nen*nen))
    col = np.zeros(int(Ne*nen*nen))

    counter = 0
    for e in range(Ne): # loop index in [0,Ne-1]
        for idx_a in range(nen): # loop index in [0,nen-1]
            for idx_b in range(nen): # loop index in [0,nen-1]
                row[counter] = ien[e,idx_a]
                col[counter] = ien[e,idx_b]
                counter = counter + 1

    assert(len(row)==len(col))

    return row, col

def get_csr_case2():
    # return car info: case 2
    Nex1 = 3 # num. elements in x1 dir
    Nex2 = 2 # num. elements in x2 dir
    Ne = int(Nex1*Nex2) # num. elements in total: 2D structured quad mesh
    Nn = int((Nex1+1)*(Nex2+1)) # num. nodes/vertices

    nen = 4 # a linear quad mesh
    ien = np.zeros([Ne,nen])

    for e in range(Ne): # loop index in [0,Ne-1]
        ien[e,0] = int(e/Nex1)+e
        ien[e,1] = int(e/Nex1)+e+1
        ien[e,2] = int(e/Nex1)+e+1+Nex1+1
        ien[e,3] = int(e/Nex1)+e+Nex1+1

    row, col = get_row_col(ien)

    nrows = int(Nn)
    ncols = nrows  
    assert(np.max(row)<=nrows-1)
    assert(np.max(col)<=ncols-1)

    # data with ones (not zeros)
    data = np.ones(len(row))

    return nrows, ncols, row, col, data

def get_csr_case3():
    # return car info: case 3
    Ne = 9 # num. elements in a linear tri mesh
    Nn = 9 # num. nodes/vertices

    nen = 3 # a linear tri mesh
    ien = np.zeros([Ne,nen])
  
    e = 0
    ien[e,0] = 0
    ien[e,1] = 1
    ien[e,2] = 3
  
    e = 1
    ien[e,0] = 1
    ien[e,1] = 4
    ien[e,2] = 3

    e = 2
    ien[e,0] = 1
    ien[e,1] = 2
    ien[e,2] = 4

    e = 3
    ien[e,0] = 0
    ien[e,1] = 3
    ien[e,2] = 5

    e = 4
    ien[e,0] = 3
    ien[e,1] = 7
    ien[e,2] = 5

    e = 5
    ien[e,0] = 3
    ien[e,1] = 4
    ien[e,2] = 7

    e = 6
    ien[e,0] = 2
    ien[e,1] = 8
    ien[e,2] = 4

    e = 7
    ien[e,0] = 5
    ien[e,1] = 7
    ien[e,2] = 6

    e = 8
    ien[e,0] = 4
    ien[e,1] = 8
    ien[e,2] = 7

    row, col = get_row_col(ien)

    nrows = int(Nn) # same as Nn
    ncols = nrows
    assert(np.max(row)<=nrows-1)
    assert(np.max(col)<=ncols-1)

    # data with ones (not zeros)
    data = np.ones(len(row))

    return nrows, ncols, row, col, data

def plt_csr_mat():

    nrows, ncols, row, col, data = get_csr_case2()

    mat = csr_matrix((data, (row, col)), shape = (nrows, ncols))
    plt.spy(mat)
    plt.savefig('sparse_mat_csr_v1_plot1.pdf')
    plt.show()

plt_csr_mat()