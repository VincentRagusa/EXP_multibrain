from pathlib import Path
import csv
import copy
import matplotlib.pyplot as plt
# -----------------------------------
# "parameters"

path_root = "../MABE/"
trait = "score_AVE"

# -----------------------------------
# helper functions

def C_ij(i, j):
    global current_generation
    return int(i in current_generation[j][1].split(","))

def X_i(i):
    global previous_generation
    return float(previous_generation[i][0])

def X_j(j):
    global current_generation
    return float(current_generation[j][0])

def cov(D1, D2, m_1, m_2):
    cov = 0
    for key in D1:
        cov += (D1[key] - m_1) * (D2[key] - m_2)
    return cov / len(D1)

# -----------------------------------
files = list(Path(path_root).glob("snapshot_data_*.csv"))     # load all data files
files = [str(f) for f in files]                               # convert path objects to strings
files.sort(key=lambda f: int(f.split(".")[-2].split("_")[-1])) # Sort into correct order.

column_index_from_key = {}
previous_generation = {}
current_generation = {}
current_generation["-1"] = ['0', ''] # add the projenator

LHS = []
RHS_1 = []
RHS_2 = []
RHS_3 = []
RHS = []

for fi, fp in enumerate(files):
    print("Loading ", fp)

    with open(fp,'r') as fileReader:

        if fi == 0: # if first file
            header = fileReader.readline().strip().split(",")
            for ki, key in enumerate(header):
                column_index_from_key[key] = ki # store column IDs alongside column keys
            print("Detected ", len(column_index_from_key), " data columns!")
        else: #all other files
            trash = fileReader.readline() # remove header
        
        csvReader = csv.reader(fileReader) # use csv to parse the rest of the file correctly

        previous_generation = copy.deepcopy(current_generation)
        current_generation = {}
        for line in csvReader:
            current_generation[line[column_index_from_key["ID"]]] = [ line[column_index_from_key[trait]], line[column_index_from_key["snapshotAncestors_LIST"]] ]
    # Finished loading data (for this file...)

    # Make Ancester -> Decendant connection matrix and helper stuff
    na = len(previous_generation)
    nd = len(current_generation)
    Cad = {i : {j : C_ij(i,j) for j in current_generation} for i in previous_generation}
    Ca = {i: sum([Cad[i][j] for j in current_generation]) for i in previous_generation}
    Cd = {j: sum([Cad[i][j] for i in previous_generation]) for j in current_generation}
    C = sum(sum([Cad[i][j] for i in previous_generation]) for j in current_generation)

    # do price stuff (Kerr & Godfrey-Smith 2009)
    ave_Ca = sum([Ca[i] for i in Ca]) / len(Ca)
    ave_Cd = sum([Cd[j] for j in Cd]) / len(Cd)

    Xa = {i: X_i(i) for i in previous_generation}
    Xd = {j: X_j(j) for j in current_generation}

    ave_Xa = sum(Xa.values()) / na
    ave_Xd = sum(Xd.values()) / nd

    delta_ave_X = ave_Xd - ave_Xa #LHS of price equation!
    # print("LHS", delta_ave_X)
    LHS.append(delta_ave_X)

    cov_Ca_Xa = cov(Ca, Xa, ave_Ca, ave_Xa) # RHS part 1 !
    # print("\tPart 1", cov_Ca_Xa / (C/na))
    RHS_1.append(cov_Ca_Xa / (C/na))

    delta_Xa = {i: sum([Cad[i][j]*(Xd[j] - Xa[i]) for j in current_generation]) / Ca[i] if Ca[i] != 0 else 0 for i in previous_generation}
    ave_Ca_delta_Xa = sum([Ca[i]*delta_Xa[i] for i in previous_generation]) / na # RHS part 2 !
    # print("\tPart 2", ave_Ca_delta_Xa / (C/na))
    RHS_2.append(ave_Ca_delta_Xa / (C/na))

    cov_Cd_Xd = cov(Cd, Xd, ave_Cd, ave_Xd) # RHS part 3 !
    # print("\tPart 3", -(nd/na)*cov_Cd_Xd / (C/na))
    RHS_3.append(-(nd/na)*cov_Cd_Xd / (C/na))

    price_RHS = (cov_Ca_Xa + ave_Ca_delta_Xa - (nd/na)*cov_Cd_Xd) / (C/na)
    # print("RHS", price_RHS)
    RHS.append(price_RHS)

print("DONE! Data loading is complete.")

import matplotlib.animation as animation

fig, axs = plt.subplots(1,1)

def plot(f):
    # global ani
    # if f > len(LHS):
    #     ani.event_source.stop()
    axs.clear()
    axs.title.set_text(f)
    axs.scatter(RHS_1[:f], LHS[:f], label="RHS 1 cov(Ca,Xa)", alpha=0.25)
    axs.scatter(RHS_2[:f], LHS[:f], label="RHS 2 ave[Ca*dXa]", alpha=0.25, marker="s")
    axs.scatter(RHS_3[:f], LHS[:f], label="RHS 3 -cov(Cd, Xd)", alpha=0.25, marker="^")
    axs.scatter(RHS[:f], LHS[:f], color="black", alpha=0.25)
    axs.legend()

    # axs.plot([RHS_1[:f], RHS_2[:f]],[LHS[:f], LHS[:f]])

# ani = animation.FuncAnimation(fig, plot, interval=0)
plot(len(LHS))

plt.tight_layout()
plt.show()