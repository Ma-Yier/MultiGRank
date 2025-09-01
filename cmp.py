import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sco


def geneRank(filename):
    geneList = []
    with open(filename, "r") as f:
        for line in f:
            geneList.append(line.split(":")[0])
    return geneList

def cmpList(l1, l2):
    n = min([len(l1), len(l2)])
    same = []
    for i in range(n):
        same.append(len(set(l1[:i+1]).intersection(set(l2[:i+1]))))
    return same

def cmpList2(l1, l2, dict1):
    n = min([len(l1), len(l2)])
    same = []
    ess = 0
    for i in range(n):
        same.append(len(set(l1[:i+1]).intersection(set(l2[:i+1]))))
        if l2[i] in dict1:
            ess += 1
        same[i] += ess
    return same

def geneEssen(modelname):
    if os.path.exists(f"{modelname}/essenlist.pkl"):
        with open(f"{modelname}/essenlist.pkl", "rb") as f:
            essenList = pickle.load(f)
            return essenList
    readm = sco.loadmat(f'{modelname}.mat', struct_as_record=False)
    matmodel = readm[modelname][0][0]
    essenList = dict()
    for  i in range(matmodel.genes.shape[0]):
        essenList[str(matmodel.genes[i][0][0])] = 1
    with open(f"{modelname}/countAll.txt", "r") as f:
        for line in f:
            gene = line.split(":")[0]
            del essenList[gene]
    print(f"essen genes: {len(essenList)}")
    with open(f"{modelname}/essenlist.pkl", "wb") as f:
        pickle.dump(essenList, f)
    return essenList

if __name__ == "__main__":
    model = "iMM904"
    essenList = geneEssen(model)
    mlen = len(model)
    filelist = os.listdir()
    genelist = dict()
    for file in filelist:
        if len(file) >= mlen + 1 and file[:mlen+1] == model + "_":
            genelist[file.split(".")[0].split("_")[-1]] = geneRank(file)
        if file == model:
            genelist["all"] = geneRank(file+"/countAll.txt")
            genelist["met"] = geneRank(file+"/countMet.txt")
    all = dict()
    met = dict()
    for key, value in genelist.items():
        if key == "all" or key == "met":
            continue
        all[key] = cmpList2(genelist["all"], value, essenList)
        met[key] = cmpList2(genelist["met"], value, essenList)
    if not os.path.exists(f"_{model}_cmp2.pkl"):
        with open(f"_{model}_cmp2.pkl", "wb") as f:
            pickle.dump({"all":all, "met":met}, f)

    plt.figure(figsize=(10, 6))
    for key, value in met.items():
        data = [value[i]/(i+1)*100 for i in range(len(value))]
        if key == "frequency":
            key = "multiplicity"
        elif key == "propotion":
            key = "frequency"
        elif key == "score":
            key = "logic"
        elif key == "naive":
            key = "degree"
        elif key == "reverse":
            key = "revdegree"
        plt.plot(data, label=key)

    #plt.title("Similarity Ratio with MetNetComp in Multiplicity Counts on " + model)
    plt.title("Similarity Ratio with MetNetComp in Frequency Counts on " + model)
    plt.xlabel("Top k Genes")
    plt.ylabel("Similarity Ratio (%)")
    plt.legend()
    plt.show()
