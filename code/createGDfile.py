import numpy as np
import matplotlib.pyplot as plt

L = 100
N = 1000

f = open('normGD.txt', "w")
f.close()

f = open('normGD.txt', "a")

Yx = {}
Yy = {}

for j in range(0, N):
    result = 0
    Yx[0], Yy[0] = np.loadtxt('guessResults/guessResults' + str(j) + '.txt', delimiter=' ', unpack=True, max_rows=L)
    Yx[1], Yy[1] = np.loadtxt('guessResults/guessResults' + str(j+1) + '.txt', delimiter=' ', unpack=True, max_rows=L)
    for i in range(L):
        result += (Yx[0][i] - Yx[1][i])**2 + (Yy[0][i] - Yy[1][i])**2
    result = np.sqrt(result)
    f.write(str(result) + '\n')

f.close()
