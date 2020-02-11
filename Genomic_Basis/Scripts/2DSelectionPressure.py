import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as LA

import matplotlib.animation as animation

fig, axs = plt.subplots(1,2)

def Fa(x):
    return x

def Fb(x):
    return -0.1*x

def C(x, la, lb):
    ca = 0
    cb = 0
    while x >= la+lb:
        ca += 1
        cb += 1
        x -= la+lb
    if x >= la:
        ca += 1
        x -= la
    return ca, cb, x

def F(x):
    la = 2
    lb = 0.5
    ca, cb, d = C(x, la, lb)
    if ca == cb:
        return Fa(ca*la + d) + Fb(cb*lb)
    else:
        return Fa(ca*la) + Fb(cb*lb + d)



def animate(f):
    names_list = ["", "A::", "B::", "C::", "D::", "E::", "F::","G::","H::","I::","J::","K::","L::","M::","N::","O::","P::","Q::","R::","S::","T::","U::","V::","W::","X::","Y::","Z::"]

    axs[0].clear()
    axs[1].clear()
    datamap = {}
    with open("../MABE/pop.csv", 'r') as inputfile:
        X = {}
        
        found_names = set()
        for i, colname in enumerate(inputfile.readline().strip().split(",")):
            for nameprefix in names_list:
                if colname == nameprefix+"mean_AVE":
                    datamap[colname] = i
                    found_names.add(nameprefix)

        for nameprefix in found_names:
            X[nameprefix] = []

        for line in inputfile:
            for nameprefix in found_names:
                X[nameprefix].append(float(line.strip().split(",")[datamap[nameprefix+"mean_AVE"]]))

    XMIN = min([min(X[key]) for key in X])
    XMAX = max([max(X[key]) for key in X])


    # for nameprefix in sorted(found_names):
    #     if nameprefix == "": continue         
    #     axs[0].plot(X[nameprefix],Y[nameprefix], label="."+nameprefix, alpha = 0.75)

    for nameprefix in sorted(found_names):
        if nameprefix == "": continue
        # if nameprefix != "A::":
            # axs[0].plot(X[nameprefix][-1],X[nameprefix][-1]/4, label="."+nameprefix, marker="^")
        # else:
        axs[0].plot(X[nameprefix][-1],F(X[nameprefix][-1]), label="."+nameprefix, marker="o")
        
    for nameprefix in sorted(found_names):
        if nameprefix == "": continue
        # if nameprefix != "A::":
        #     xx = list(range(2,len(X[nameprefix])))
        #     yy = [X[nameprefix][i-1]/4-X[nameprefix][i-2]/4 for i in xx]
        #     axs[1].plot(xx, yy, label="."+nameprefix)
        # else:
        xx = list(range(2,len(X[nameprefix])))
        yy = [F(X[nameprefix][i-1])-F(X[nameprefix][i-2]) for i in xx]
        axs[1].plot(xx, yy, label="."+nameprefix)
    axs[1].axhline(color="black")

    x = np.arange(XMIN*0.9,XMAX*1.1, 0.01)
    y = [F(xx) for xx in x]
    axs[0].plot(x,y, label="fitness function")
    # axs[0].plot(x,np.divide(x,4), label="other function")

    axs[0].legend()

ani = animation.FuncAnimation(fig, animate, interval=1)
plt.tight_layout()
plt.show()