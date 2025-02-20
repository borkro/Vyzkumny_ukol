import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11

L = 30

X = np.linspace(0, L-1, L)
Y = np.loadtxt('costFunction.txt', unpack=True, max_rows=L)

plt.plot(X, Y, 'ro', linewidth=1, linestyle='solid')
# plt.title('Účelová funkce J')
plt.xlabel('Iterace')
plt.ylabel('J')
# plt.yscale('log')
# plt.axis([0.0, L-1, 0.0, 3700.0])
plt.tight_layout()
plt.grid()
# plt.savefig('plotAdjoint.png')
plt.show()
