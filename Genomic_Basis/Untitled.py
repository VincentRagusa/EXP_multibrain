#!/usr/bin/env python
# coding: utf-8

# In[8]:


import matplotlib.pyplot as plt
plt.figure(figsize=(20,20))

datamap = {}
with open("MABE/pop.csv", 'r') as inputfile:
    Xx = {}
    Yy = {}
    for nameprefix in ["", "A::", "B::", "C::", "D::", "E::", "F::"]:
        Xx[nameprefix] = []
        Yy[nameprefix] = []
        
    for i, colname in enumerate(inputfile.readline().strip().split(",")):
        for nameprefix in ["", "A::", "B::", "C::", "D::", "E::", "F::"]:
            if colname in [nameprefix+"variance_AVE", nameprefix+"mean_AVE"]:
                datamap[colname] = i

    for j, line in enumerate(inputfile):
        if j == 0: continue
        for nameprefix in ["", "A::", "B::", "C::", "D::", "E::", "F::"]:
            Xx[nameprefix].append(float(line.strip().split(",")[datamap[nameprefix+"variance_AVE"]]))
            Yy[nameprefix].append(float(line.strip().split(",")[datamap[nameprefix+"mean_AVE"]]))
#         if j > 1000: break
            
for nameprefix in ["", "A::", "B::", "C::", "D::", "E::", "F::"]:         
    plt.plot(Xx[nameprefix],Yy[nameprefix], label="."+nameprefix)
plt.legend()
# plt.show()


# In[19]:


import math
import numpy as np
x = list(range(0,60))
y = list(range(0,60))
def fn(x,y):
    return np.sqrt(x)* np.sin(np.sin(y)*x) + 8*y

#4*np.sin(4*x)*
# In[ ]:





# In[21]:


# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


fig = plt.figure(figsize=(20,20))
ax = fig.gca(projection='3d')

# Make data.
X = np.arange(20, 38, 0.05)
Y = np.arange(20, 38, 0.05)
X, Y = np.meshgrid(X, Y)
R = fn(X,Y)
Z = R

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False, rcount=200,ccount=200)

# Customize the z axis.
# ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
# plt.show()

# fig = plt.figure(figsize=(20,20))
# ax = fig.add_subplot(111, projection='3d')
for nameprefix in ["", "A::", "B::", "C::", "D::", "E::", "F::"]:
    print(nameprefix)
    ax.plot(Xx[nameprefix],Yy[nameprefix], zs = [fn(p[0],p[1]) for p in zip(Xx[nameprefix],Yy[nameprefix]) ], label="."+nameprefix, linewidth = 5)

plt.xlim(32,38)
plt.ylim(32,38)
plt.show()
# rotate the axes and update
# for angle in range(0, 360):
    # ax.view_init(30, 45)
    # plt.draw()
    # plt.pause(.001)


# In[ ]:




