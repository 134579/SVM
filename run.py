import quadprog
from cvxopt import matrix
import numpy as np
import matplotlib.pyplot as plt
import cvxopt
import svm

def drawline(w, b):
	x = np.linspace(-1, 11, 100)
	x = np.reshape(x, (100, 1))
	y = (-b - w[0] * x) / w[1]
	return [x, y]


# number of points
N = 50
cvxopt.setseed(0)
X = cvxopt.uniform(2, N, 0, 10)
# create a line crossing point (5,5)
w = cvxopt.uniform(2, 1, 0, 1)
b = 0 - w.T * matrix([5, 5])
label = []
for i in range(N):
	if (w.T * X[:, i] + b)[0, 0] > 0:
		label.append(1)
	else:
		label.append(-1) 
plt.figure('1')
plt.scatter(X[0, :], X[1, :], c=label, cmap='cool')
[x, y] = drawline(w, b)
plt.plot(x, y, label='original line')



[w, b] = svm.train(X, label);

[x, y] = drawline(w, b)
plt.plot(x, y, label='fitted line')
plt.legend()
plt.show()






