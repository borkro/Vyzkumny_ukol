from tkinter.ttk import Separator
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}\usepackage{float}\usepackage{textcomp}\usepackage{mathtools}\usepackage{amsmath}\usepackage{amsthm}\usepackage{amssymb}\usepackage{graphicx}\usepackage{setspace}\renewcommand{\vec}[1]{\boldsymbol{#1}}\usepackage[varg]{txfonts}'

N = 1000

X = np.linspace(0, N - 1, N)

norm = np.loadtxt('normGD.txt', delimiter=' ', unpack=True, max_rows=N)


plt.plot(X, norm, 'r-', linewidth=1, linestyle='solid')
plt.xlabel(r'$k$')
plt.ylabel(r'$l_2$ norm')
plt.grid(axis='y')
plt.title(r'Norm of GD step')
#plt.legend()

plt.tight_layout()
plt.show()
plt.savefig('normPlot.eps', bbox_inches='tight')