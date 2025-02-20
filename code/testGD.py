import numpy as np
import matplotlib.pyplot as plt

Yx = {}
Yy = {}
L = 50
Yx[0], Yy[0] = np.loadtxt('guessResults/guessResults0.txt', delimiter=' ', unpack=True, max_rows=L)
Yx[1], Yy[1] = np.loadtxt('guessResults/guessResults1000.txt', delimiter=' ', unpack=True, max_rows=L)

result = 0
for i in range(L):
    result += (Yx[0][i] - Yx[1][i])**2 + (Yy[0][i] - Yy[1][i])**2
result = np.sqrt(result)
print(result)