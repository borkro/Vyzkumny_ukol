from tkinter.ttk import Separator
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}\usepackage{float}\usepackage{textcomp}\usepackage{mathtools}\usepackage{amsmath}\usepackage{amsthm}\usepackage{amssymb}\usepackage{graphicx}\usepackage{setspace}\renewcommand{\vec}[1]{\boldsymbol{#1}}\usepackage[varg]{txfonts}'

N = 100

X, Y = np.loadtxt('linePlot.txt', delimiter=' ', unpack=True, max_rows=N)


plt.plot(X, Y, 'r-', linewidth=1, linestyle='solid')
plt.xlabel(r'Step')
plt.ylabel(r'Cost function')
plt.grid(axis='y')
plt.title(r'Cost function on line')
#plt.legend()

plt.tight_layout()
plt.savefig('plotLine.eps', bbox_inches='tight')