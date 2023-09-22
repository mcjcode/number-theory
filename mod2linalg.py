import numpy as np


def solve_linear_equations_mod_2(A, b):
    """
    Solves a system of linear equations over the field of integers
    mod 2.

    Args:
      A: A NumPy array of the coefficients of the system of equations.
      b: A NumPy array of the values on the right-hand side of the equations.

    Returns:

    A NumPy array of the solutions to the system of equations, or None
    if there is no solution.

    """

    origA = A.copy()
    origb = b.copy()
    
    b = b[:,np.newaxis]
    A = np.hstack((A, b))
    # iterate over the columns of A
    for i in range(A.shape[0]):
        # Find the pivot row
        for j in range(i, A.shape[0]):
            if A[j, i]:
                break
        # Swap the current row with the pivot row.
        if i != j:
            A[[i, j]] = A[[j, i]]
        # Subtract the current row from all the other rows
        for j in range(A.shape[0]):
            if i != j:
                if A[j,i]:
                    A[j,:] = A[j,:] ^ A[i,:]    
    # The solution vector is now stored in the last column of A.
    if (origA.dot(A[:,-1])%2 == origb).all():
        return A[:, -1]
    else:
        raise ValueError('Solve Failed')
    return A[:, -1]


def null_space_mod_2(A):
    """
    Return a basis for the null space
    of the matrix A.
    """
    #
    # first row reduce the matrix A
    #
    A = A.copy()
    
    i0 = 0
    for j in range(A.shape[1]):
        for i in range(i0, A.shape[0]):
            if A[i, j]:
                if i != i0:
                    A[[i0, i]] = A[[i, i0]]
                for i1 in range(A.shape[0]):
                    if i1 != i0:
                        if A[i1, j]:
                            A[i1,:] = A[i1,:] ^ A[i0,:]
                i0 += 1
                break
        if i0==A.shape[0]:
            break

    pivots = []
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if A[i,j]:
                pivots.append((i,j))
                break    
    pivot_cols = [j for (_, j) in pivots]
    non_pivot_cols = [j for j in range(A.shape[1]) if j not in pivot_cols]
    
    basis = []
    for j in non_pivot_cols:
        v = np.zeros((A.shape[1],), dtype=np.int8)
        v[j] = 1
        for (i0, j0) in reversed(pivots):
            v[j0] = v[j0+1:].dot(A[i0,j0+1:]) % 2
        basis.append(v)
    if len(basis)==0:
        return np.zeros((A.shape[1], 0), dtype=np.int8)
    else:
        return np.array(basis).transpose()


def det_mod_2(A):
    for i in range(A.shape[0]):
        # Find the pivot row
        for j in range(i, A.shape[0]):
            if A[j, i]:
                break
        if i != j:
            A[[i,j]] = A[[j,i]]
        for j in range(A.shape[0]):
            if j>i:
                if A[j,i]:
                    A[j,:] = A[j,:] ^ A[i,:]
    return np.prod(A.diagonal())
