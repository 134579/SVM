# X : D x N  sample set, cvxopt.matrix object
# label : list
import numpy as np
from cvxopt import matrix  
import quadprog
def train(X, label):
	[D, N] = X.size
	H = matrix(np.identity(D + 1))
	H[D, D] = 0;
	g = matrix(np.zeros((D + 1, 1)))
	A2 = matrix(np.zeros((N, D + 1)))
	for i in range(N):
		for j in range(D):
			A2[i, j] = label[i] * X[j, i]
		A2[i, D] = label[i] * 1
			
	A2 = A2.T
	b2 = matrix(np.ones((N, 1)))
	wb = quadprog.quadprog(H, g, A2=A2, b2=b2)
	w = wb[0:D]
	b = wb[D:D + 1]
	return [w, b]
