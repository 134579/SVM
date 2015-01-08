import numpy  as np
from cvxopt import matrix, solvers
import cvxopt
from __builtin__ import len, min
from copy import copy

MIN = 1e-5

def allzero(A):
	[m, n] = A.size
	for i in range(m):
		for j in range(n):
			if abs(A[i, j]) >= MIN:
				return False
	return True

def allgezero(A):
	[m, n] = A.size
	for i in range(m):
		for j in range(n):
			if A[i, j] < -MIN:
				return False
	return True
	
# A^T*x = b
# returns [] : I	
def getI(A, b, x):
	assert A.size[0] == x.size[0], 'Usage: get i where A[:,i].T*x = b_i'
	m = A.size[1]
	r = []
	for i in range(m):
		if allzero(A[:, i].T * x - b[i, 0]):
			r.append(i)
	return r
			


# min Q(x) = 1/2*x^T*H*x + g^T*x 
# st. A1^T*x = b1
#     A2^T*x >=b2
# !! all of the matrixes are cvxopt.matrix
#   not np.ndarray
def quadprog(H, g, A1=0, b1=0, A2=0, b2=0):
	if isinstance(A2, int):
		quadprogWithEquation(H, g, A1, b1)
	elif isinstance(A1, int): 
		# for svm
		assert H.size[0] == H.size[0], 'H must be a symetric matrix'
		assert H.size[0] == A2.size[0]  
		H = 1.00 * H
		g = 1.00 * g
		A2 = 1.00 * A2
		b2 = 1.00 * b2
		n = H.size[0]  # dimension of x
		m = A2.size[1]  # # of inequalities
		
		
		N = 1000
		x = range(N)
		S = range(N)
		d = range(N)
		alpha = range(N)
		Lambda = range(N)
		
		# step 1, calculate first feasible point
		c = A2[:, 0]
		tmp = solvers.lp(c, -1 * A2.T, -1 * b2)
		x[0] = tmp['x']
		# x.append(tmp['x'])  # x[0]
		I = getI(A2, b2, x[0])
		S[0] = I
		# S.append(I)  # S[0]
		
		
		k = 0
		while k < N - 1:
			# step 2
			# min 1/2*d.T*H*d + (g+H.T*xk).T*d
			# s.t A[:,i].T*d = 0 where i \in I
			zeros = matrix(np.zeros((len(S[k]), 1))) * 1.0
			d_and_lambda = quadprogWithEquation(H, g + H.T * x[k], A2[:, S[k]], zeros)
			d[k] = d_and_lambda[0:n, 0]
			Lambda[k] = d_and_lambda[n:, 0]
			# d.append(d_and_lambda[0:n, 0])
			# Lambda.append(d_and_lambda[n:, 0])
									
			if not allzero(d[k]):  # goto step 3
				# step 3, calculate alphak
				alpha_cand = [1]
				for i in range(m):
					if S[k].count(i) == 0:
						if (A2[:, i].T * d[k][0, 0])[0, 0] < 0:
							alpha_cand.append((b2[i] - A2[:, i].T * x[k])[0, 0] / (A2[:, i].T * d[k])[0, 0])
				# alpha.append(min(alpha_cand))
				alpha[k] = min(alpha_cand)
				# x.append(x[k] + alpha[k] * d[k])
				x[k + 1] = x[k] + alpha[k] * d[k]
				if alpha[k] == 1:
					# step 4
					# S.append(copy(S[k]))
					S[k + 1] = copy(S[k])
					k = k + 1
					continue
				else:
					# add a constraint
					I = getI(A2, b2	, x[k] + alpha[k] * d[k])
					for xx in I:
						if xx not in S[k]:
							S[k].append(xx)
					# step 4
					S[k + 1] = copy(S[k])
					# S.append(copy(S[k]))
					k = k + 1
					continue
			else:
				if allgezero(Lambda[k]):
					# exit
					return x[k]
				else:
					ik = 0
					for xx in range(Lambda[k].size[0]):
						if Lambda[k][xx] < Lambda[k][ik]:
							ik = xx
					S[k].remove(S[k][ik])
					# x.append(copy(x[k]))
					x[k + 1] = copy(x[k])
					# step 4
					# S.append(copy(S[k]))
					S[k + 1] = copy(S[k])
					k = k + 1
					continue

		return x[N - 1]

# min Q(x) = 1/2*x^T*H*x + g^T*x 
# st. A^T*x = b
# Return : [x,lambda].T 
# !! only work for cvxopt.atrix class

def quadprogWithEquation(H, g, A, b):
	quadprogWithEquation.xxxxx = quadprogWithEquation.xxxxx + 1
	print quadprogWithEquation.xxxxx
	if A.size[0] != 0:
		assert H.size[0] == H.size[0], 'H must be a symetric matrix'
		assert H.size[0] == g.size[0]
		assert H.size[0] == A.size[0]
		assert b.size[0] == A.size[1]  
		# just solve:
		# g + H*x = A*lambda
		# A^T*x = b
		m = A.size[1]
		
		T = -1.00 * matrix(np.vstack((g, b)))
		
		S = matrix(np.vstack((np.hstack((H, -1 * A)), np.hstack((-1 * A.T, np.zeros((m, m)))))))
	
		tmp = matrix(T)
		cvxopt.lapack.gesv(S, tmp)
		return tmp
	else:
		# just solve:
		# g + H*x = 0
		assert H.size[0] == H.size[0] , 'H must be a symetric matrix'
		assert H.size[0] == g.size[0]
		S = 1.00 * H
		T = -1.00 * g
		tmp = matrix(T)
		cvxopt.lapack.gesv(S, tmp)
		return tmp


quadprogWithEquation.xxxxx = 0
