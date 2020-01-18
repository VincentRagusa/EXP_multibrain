import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as LA

import matplotlib.animation as animation

fig, axs = plt.subplots(1,2)

selection_pressure_history = {}
sph_init = True

def fn(x,y):
    #return np.sin(y)*np.sin(x-y) + y
    # return 4*np.sin(x-y/4) + 2*y
    # return np.sqrt(x) * np.sin(np.sin(y)*x) + 2*y
    # return np.sqrt(x)*np.sin(8*x)* np.sin(np.sin(y)*x) + 2*y
    # return np.sin(y) + y * 0*x
    return 0*x + y + np.sqrt(145.0)*np.cos(y)*np.sqrt(1.0/(1.0+144.0*np.power(np.cos(y),2)))
    # return x*y



def grad_fn(x,y):
    # below for np.sqrt(x)*np.sin(8*x)* np.sin(np.sin(y)*x) + 2*y
    # https://www.wolframalpha.com/input/?i=partial+derivative+sqrt%28x%29*sin%288*x%29*sin%28sin%28y%29*x%29+%2B+2*y
    # ---
    # partial_x = ((2*x*np.sin(8*x)*np.sin(y)*np.cos(x*np.sin(y))) + (np.sin(x*np.sin(y))*np.sin(8*x)) + (16*x*np.cos(8*x))) / (2*np.sqrt(x))
    # partial_y = (np.power(x, 3.0/2.0) * np.sin(8*x) * np.cos(y) * np.cos(x*np.sin(y)) ) + 2
    # ---
    # below for np.sin(x) + x
    # https://www.wolframalpha.com/input/?i=partial+derviatives+of+sin%28x%29+%2B+x
    # ---
    # partial_x = 0
    # partial_y = np.cos(y) + 1
    # ---
    # below for 0*x + y + np.sqrt(145.0)*np.cos(y)*np.sqrt(1.0/(1.0+144.0*np.power(np.cos(y),2)))
    # https://www.wolframalpha.com/input/?i=partial+derivatives+sqrt%28%281%2B12%5E2%29%2F%281%2B%2812%5E2%29*cos%5E2%28v%29%29%29*cos%28v%29+%2B+v
    # https://www.wolframalpha.com/input/?i=derivative+x+%2B+sqrt%28145%29+cos%28x%29+sqrt%281%2F%281+%2B+144+cos%5E2%28x%29%29%29
    # ---
    partial_x = 0
    partial_y = (144.0*np.sqrt(145.0)*np.sin(y)*np.power(1.0/(144.0*np.power(np.cos(y),2)+1), 3.0/2.0)*np.power(np.cos(y),2)) - (np.sqrt(145.0)*np.sin(y)*np.sqrt(1/(144.0*np.power(np.cos(y),2)+1.0))) + 1.0
    # ---
    
    print(y, partial_y)
    return (partial_x, partial_y)


def normalize_vector(vector, p=2):
    length = LA.norm(vector, p)
    return list(np.divide(vector, length))


def animate(f):
    global sph_init, selection_pressure_history
    names_list = ["", "A::", "B::", "C::", "D::", "E::", "F::","G::","H::","I::","J::","K::","L::","M::","N::","O::","P::","Q::","R::","S::","T::","U::","V::","W::","X::","Y::","Z::"]

    axs[0].clear()
    axs[1].clear()
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
            if sph_init:
                selection_pressure_history[nameprefix] = []
        sph_init = False

        for j, line in enumerate(inputfile):
            # if j % 2 == 0: continue
            for nameprefix in found_names:
                X[nameprefix].append(float(line.strip().split(",")[datamap[nameprefix+"variance_AVE"]]))
                Y[nameprefix].append(float(line.strip().split(",")[datamap[nameprefix+"mean_AVE"]]))

    YMIN = min([min(Y[key]) for key in Y])
    XMIN = min([min(X[key]) for key in X])

    YMAX = max([max(Y[key]) for key in Y])
    XMAX = max([max(X[key]) for key in X])


    for nameprefix in sorted(found_names):
        if nameprefix == "": continue         
        axs[0].plot(X[nameprefix],Y[nameprefix], label="."+nameprefix, alpha = 0.75)

    for nameprefix in sorted(found_names):
        if nameprefix == "": continue
        axs[0].plot(X[nameprefix][-1],Y[nameprefix][-1], color='white',marker="o")
        
        if nameprefix != "A::": #TODO special thing
            gradient = [0,1/2]
        else:
            gradient = grad_fn(X[nameprefix][-2],Y[nameprefix][-2])

        n_gradient = normalize_vector(gradient)
        v_gradient = np.multiply(n_gradient, min(YMAX-YMIN, XMAX-XMIN)/25)
        axs[0].plot([X[nameprefix][-1],X[nameprefix][-1]+v_gradient[0]],[Y[nameprefix][-1],Y[nameprefix][-1]+v_gradient[1]], color="red")

        history = [X[nameprefix][-1]-X[nameprefix][-2],Y[nameprefix][-1]-Y[nameprefix][-2]]
        n_history = normalize_vector(history)
        v_history = np.multiply(n_history, min(YMAX-YMIN, XMAX-XMIN)/25)
        axs[0].plot([X[nameprefix][-1],X[nameprefix][-1]+v_history[0]],[Y[nameprefix][-1],Y[nameprefix][-1]+v_history[1]], color="blue")

        selection_pressure_history[nameprefix].append(np.dot(n_gradient, n_history)*LA.norm(history, 2))#abs(fn(X[nameprefix][-1],Y[nameprefix][-1])-fn(X[nameprefix][-2],Y[nameprefix][-2]))) #

    for nameprefix in sorted(found_names):
        if nameprefix == "": continue 
        axs[1].plot(range(len(selection_pressure_history[nameprefix])), selection_pressure_history[nameprefix], label="."+nameprefix)
        # axs[1].plot(range(len(selection_pressure_history[nameprefix][-50:])), selection_pressure_history[nameprefix][-50:], label="."+nameprefix)
    axs[1].axhline(color="black")

    x = np.arange(XMIN*0.95,XMAX*1.05, XMAX/1000)
    y = np.arange(YMIN*0.95,YMAX*1.05, YMAX/1000)
    xx, yy = np.meshgrid(x, y, sparse=True)
    z = fn(xx,yy)
    axs[0].contourf(x,y,z, alpha = 0.25,levels=15)

    axs[0].legend()
    # axs[1].legend()
    # plt.show()

ani = animation.FuncAnimation(fig, animate, interval=1)
plt.tight_layout()
plt.show()