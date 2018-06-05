'''
==============
3D scatterplot
==============
'''

import importlib
importlib.import_module('mpl_toolkits').__path__
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xs = np.array(pd.read_csv('dimX.csv', sep=',',header=None))
ys = np.array(pd.read_csv('dimY.csv', sep=',',header=None))
zs = np.array(pd.read_csv('dimZ.csv', sep=',',header=None))
c = np.array(pd.read_csv('color.csv', sep=',',header=None)).astype(int)
co = []
for item in c:
    co.append(item[0])

ax.scatter(xs, ys, zs, c=co, marker='o', s = 10)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
