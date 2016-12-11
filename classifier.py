import numpy as np
from sklearn import svm
from sklearn.cross_validation import KFold
from random import uniform
from sklearn.cross_validation import cross_val_score

## data is N patient samples x 251 pathweights (corresponding to each unique src-target pair)
## target is N patient samples x 1 (takes on value +1 if lung, -1 if colon)

## make some dummy data
data = np.zeros((60, 251))
for i in xrange(60):
	for j in xrange(251):
		data[i][j] = uniform(0, 100)
target = np.zeros(60)
for i in xrange(60):
	if i < 30:
		target[i] = 1
	else:
		target[i] = -1

## try to do classification with 5-fold classification
kf = KFold(n=60, n_folds=5)
model = svm.SVC(kernel='linear', C=1)
results = cross_val_score(model, data, target, cv=kf)
print np.mean(results)
print np.std(results)