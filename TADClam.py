import numpy as np
 

A = np.loadtxt("raw_chr18_300_500_30kb.tsv")
U, S, Vh = np.linalg.svd(A, full_matrices=False)
arr = np.zeros((50,50))
for i in range(len(S)):
    arr[i,i]=S[i]
np.savetxt("aftersvd.txt",U.dot(arr).dot(Vh))
 
