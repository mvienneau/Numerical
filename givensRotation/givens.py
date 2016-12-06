#! /usr/bin/env python3

import numpy


def givens_Q_times(V, x):
    """Apply orthogonal matrix represented as list of givens rotations"""
    y = x.copy()
    for i in reversed(range(len(V))):
        y[i:] += numpy.dot(V[i].T, y[i:])
    return y

def rot_matrix(t):
    return numpy.array([[numpy.cos(t), numpy.sin(t)], [-numpy.sin(t), numpy.cos(t)]])

def givens_rot(theta, v, x):
    rotation = rot_matrix(theta)
    return rotation

def givens_rot2(theta, v, x):
    retMatrix = x.copy()
    m,n = x.shape
    U1 = numpy.zeros(m)
    U1[0] = 1
    U2 = numpy.zeros(m)
    for i in range(0,m):
        if (i == 0):
            U2[0] = 0
        else:
            U2[i] = v[i-1]
    U = numpy.array([U1.T, U2.T]).T
    G = rot_matrix(theta)
    rot = (G.dot(U.T) - U.T)
    rot1 = U.dot(rot)
    """going through each column """
    for i in range(0, n):
        y = retMatrix.T[i]
        """getting the whole column of the matrix 'x' """
        y = U.dot((G.dot(U.T).dot(y) - U.T.dot(y)))
        retMatrix.T[i] = y
    return retMatrix, rot1

def givens_main2(A):
    V = []
    G = []
    m,n = A.shape
    R = A.copy()
    for i in range(0,n):
        """x is a column vector of R """
        x = R[i:, i].copy()
        norm = (numpy.linalg.norm(x[1:]))
        """v is the normalized column vector from one down (below diaganol) """
        v = x[1:]/norm
        """theta is the angle in which you rotate X to introduce all zeros
        also rotate the rest of the columns in that same way (will not introduce special zeros, however) """
        theta = numpy.arctan2(norm, x[0])
        """ add (theta,v) to the list in order to repeady apply this givens rot method to the identity matrix """
        r, rot = givens_rot2(theta, v, R[i:, i:])
        R[i:, i:] += r
        V.append(rot)
    Q = numpy.eye(m,n)
    R = numpy.triu(R[:n,:])
    for i in range(n):
        Q[:,i] = givens_Q_times(V, Q[:,i])
    return Q, R


def givens_main(A):
    V = []
    G = []
    m,n = A.shape
    for j in range(0,n):
        for i in range(m, j+1, -1):
            gmatrix = numpy.eye(m,m)
            x = A[i-2:i, j:n]
            norm = (numpy.linalg.norm(x[1:]))
            v = x[1:]/norm
            theta = numpy.arctan2(x[1,0], x[0,0])
            V.append((theta,v))
            gmatrix[i-2:i, i-2:i] = givens_rot(theta, v, x)
            G.append(gmatrix)
            A = gmatrix.dot(A)
    Q = numpy.eye(m,m)

    for i in range(len(G)):
        Q = Q.dot(G[i].T)
        #print (Q)
    return A, Q.T




A = numpy.array([
[3.,2,1],
[2,-3,4],
[5, 1,-1],
[7,4,2]])


print ("A: ", A)
Q, R = givens_main2(A)
print ("QTQ = I: ", Q.T.dot(Q))
print ("QR = A: ", Q.dot(R))

x = numpy.linspace(-1,1,10)
B = numpy.vander(x)
print ("Condition number of linspace m = 10: ", numpy.linalg.cond(B))
Q, R = givens_main2(B)
print ("Condition number of QR: ", numpy.linalg.cond(Q.dot(R)))
print ("Condition number of QTQ: ", numpy.linalg.cond(Q.T.dot(Q)))

x = numpy.linspace(-1,1,20)
B = numpy.vander(x)
print ("Condition number of linspace m = 20: ", numpy.linalg.cond(B))
Q, R = givens_main2(B)
print ("Condition number of QR: ", numpy.linalg.cond(Q.dot(R)))
print ("Condition number of QTQ: ", numpy.linalg.cond(Q.T.dot(Q)))
