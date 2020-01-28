import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as LA

import matplotlib.animation as animation

fig, axs = plt.subplots(1,1)


def animate(f):
    names_list = ["", "A::", "B::", "C::", "D::", "E::", "F::","G::","H::","I::","J::","K::","L::","M::","N::","O::","P::","Q::","R::","S::","T::","U::","V::","W::","X::","Y::","Z::"]

    axs.clear()
    datamap = {}

    with open("../MABE/pop.csv", 'r') as inputfile:
        X = {}
        
        found_names = set()
        for i, colname in enumerate(inputfile.readline().strip().split(",")):
            for nameprefix in names_list:
                if colname == nameprefix+"score_AVE":
                    datamap[colname] = i
                    found_names.add(nameprefix)

        for nameprefix in found_names:
            X[nameprefix] = []

        for line in inputfile:
            for nameprefix in found_names:
                X[nameprefix].append(float(line.strip().split(",")[datamap[nameprefix+"score_AVE"]]))

    
    for nameprefix in found_names:
        if nameprefix == "": continue
        traits = len(X) - 1 #-1 to account for the ave
        gens = len( X[nameprefix] )
        SC = []
        for g in range(gens):
            if X[""][g] != 0 and g > 0:
                del_trait = X[nameprefix][g] - X[nameprefix][g-1]
                del_org = X[""][g] - X[""][g-1]
                SC.append( del_trait / (del_org * traits) if del_org != 0 else 0)
        axs.plot(range(len(SC)), SC, label="."+nameprefix)

    POPSIZE = 100

    axs.plot([0,gens], [1/POPSIZE, 1/POPSIZE], label="1/N", color="black")
    axs.plot([0,gens], [-1/POPSIZE, -1/POPSIZE], label="-1/N", color="red")

# ani = animation.FuncAnimation(fig, animate, interval=1)
animate(0)
plt.legend()

plt.show()