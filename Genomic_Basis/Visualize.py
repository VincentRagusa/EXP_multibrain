import matplotlib.pyplot as plt
import numpy as np

import matplotlib.animation as animation

def fn(x,y):
    return np.sqrt(x)*np.sin(8*x)* np.sin(np.sin(y)*x) + 2*y
    # return x*y

fig, axs = plt.subplots(1,1)

def animate(f):
    names_list = ["", "A::", "B::", "C::", "D::", "E::", "F::","G::","H::","I::","J::","K::","L::"]
    #names_list = [""]
    axs.clear()
    datamap = {}
    with open("MABE/pop.csv", 'r') as inputfile:
        X = {}
        Y = {}
        
        found_names = set()
        for i, colname in enumerate(inputfile.readline().strip().split(",")):
            for nameprefix in names_list:
                if colname in [nameprefix+"variance_AVE", nameprefix+"mean_AVE"]:
                    datamap[colname] = i
                    found_names.add(nameprefix)

        for nameprefix in found_names:
            X[nameprefix] = []
            Y[nameprefix] = []

        for j, line in enumerate(inputfile):
            if j == 0: continue
            for nameprefix in found_names:
                X[nameprefix].append(float(line.strip().split(",")[datamap[nameprefix+"variance_AVE"]]))
                Y[nameprefix].append(float(line.strip().split(",")[datamap[nameprefix+"mean_AVE"]]))
    
    YMAX = max([max(Y[key]) for key in Y])
    XMAX = max([max(X[key]) for key in X])


    for nameprefix in sorted(found_names):         
        axs.plot(X[nameprefix],Y[nameprefix], label="."+nameprefix)

    x = np.arange(0,XMAX*1.05, XMAX/2000)
    y = np.arange(0,YMAX*1.05, YMAX/2000)
    xx, yy = np.meshgrid(x, y, sparse=True)
    z = fn(xx,yy)
    axs.contourf(x,y,z, alpha = 0.25,levels=30)

    # axs.xlim(0,XMAX*1.05)
    # axs.ylim(0,YMAX*1.05)
    axs.legend()
    # axs.tight_layout()


ani = animation.FuncAnimation(fig, animate, interval=1000)
plt.show()