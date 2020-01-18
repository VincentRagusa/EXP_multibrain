import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as LA

import matplotlib.animation as animation

def fn(x,y):
    #return np.sin(y)*np.sin(x-y) + y
    # return 4*np.sin(x-y/4) + 2*y
    # return np.sqrt(x) * np.sin(np.sin(y)*x) + 2*y
    return np.sqrt(x)*np.sin(8*x)* np.sin(np.sin(y)*x) + 2*y
    # return x*y

fig, axs = plt.subplots(1,1)

def animate(f):
    names_list = ["", "A::", "B::", "C::", "D::", "E::", "F::","G::","H::","I::","J::","K::","L::","M::","N::","O::","P::","Q::","R::","S::","T::","U::","V::","W::","X::","Y::","Z::"]
    #names_list = [""]
    axs.clear()
    datamap = {}
    with open("../MABE/pop.csv", 'r') as inputfile:
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
            # if j == 0: continue
            for nameprefix in found_names:
                X[nameprefix].append(float(line.strip().split(",")[datamap[nameprefix+"variance_AVE"]]))
                Y[nameprefix].append(float(line.strip().split(",")[datamap[nameprefix+"mean_AVE"]]))
    
    YMIN = min([min(Y[key]) for key in Y])
    XMIN = min([min(X[key]) for key in X])

    YMAX = max([max(Y[key]) for key in Y])
    XMAX = max([max(X[key]) for key in X])


    for nameprefix in sorted(found_names):         
        axs.plot(X[nameprefix],Y[nameprefix], label="."+nameprefix, alpha = 0.75)

    for nameprefix in sorted(found_names):
        axs.plot(X[nameprefix][-1],Y[nameprefix][-1], color='white',marker="o")

    x = np.arange(XMIN*0.95,XMAX*1.05, XMAX/1000)
    y = np.arange(YMIN*0.95,YMAX*1.05, YMAX/1000)
    xx, yy = np.meshgrid(x, y, sparse=True)
    z = fn(xx,yy)
    axs.contourf(x,y,z, alpha = 0.25,levels=15)

    # axs.xlim(0,XMAX*1.05)
    # axs.ylim(0,YMAX*1.05)
    axs.legend()
    # axs.tight_layout()


ani = animation.FuncAnimation(fig, animate, interval=10)
plt.show()