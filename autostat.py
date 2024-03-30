import json
import numpy as np
true=True
false=False
def extractline(line):
    line = line.strip()
    json_dict = json.loads(line)
    return json_dict

def readlog(logf):
    records= []
    with open(logf,"r") as f:
        for line in f:
            line = line.strip()
            if len(line)>0:
                # print(line)
                records.append(extractline(line))
    f.close()
    return records



def autostatgraph(perfs,basis,evalkeys):
    for k in evalkeys:
        res = np.array([perf[k] for perf in perfs])
        # res = res/basis[k]
        # res = res**2
        print("perf on {}, avg is {:.4f}, std is {:.4f}".format(k,np.mean(res),np.std(res,ddof=1)))
def autostat(perfs,evalkeys):
    for k in evalkeys:
        if k == "top_ds":
            res = np.array([perf[k][0] for perf in perfs])
        elif k == "hit":
            res = np.array([100*perf[k] for perf in perfs])
        else:
            res = np.array([perf[k] for perf in perfs])
        print("perf on {}, avg is {:.4f}, std is {:.4f}".format(k,np.mean(res),np.std(res,ddof=1)))
if __name__ == "__main__":
#     basis = {"dataname": "sbm", "testperf": true, "degree": 0.000341, "spec": 0.003485, "clustering": 0.033118, "orbit": 0.030991}
#     perfs = [
#         {"time": "2024-03-26 00:33:36.549315", "train_method": "gdpo", "innerloop": 1, "partial": false, "degree": 0.000491, "spec": 0.003534, "clustering": 0.029238, "orbit": 0.024308, "non_iso": 1.0, "uniq_frac": 1.0, "sbm_acc": 0.789062},

# {"time": "2024-03-26 00:33:36.549315", "train_method": "gdpo", "innerloop": 1, "partial": false, "degree": 0.000932, "spec": 0.004208, "clustering": 0.029141, "orbit": 0.02546, "non_iso": 1.0, "uniq_frac": 1.0, "sbm_acc": 0.792969},

# {"time": "2024-03-26 00:33:36.549315", "train_method": "gdpo", "innerloop": 1, "partial": false, "degree": 0.000477, "spec": 0.003275, "clustering": 0.029342, "orbit": 0.022998, "non_iso": 1.0, "uniq_frac": 1.0, "sbm_acc": 0.707031},

# {"time": "2024-03-26 00:33:36.549315", "train_method": "gdpo", "innerloop": 1, "partial": false, "degree": 0.000544, "spec": 0.003441, "clustering": 0.029094, "orbit": 0.021671, "non_iso": 1.0, "uniq_frac": 1.0, "sbm_acc": 0.757812},


#     ]
#     evalkeys = ["spec"]
#     autostatgraph(perfs,basis,evalkeys)
    perf = [
        {"timestamp": "03-27-19:25:09", "target": "jak2", "validity": 1.0, "uniqueness": 1.0, "novelty": 0.2981366459627329, "top_ds": [8.200000000000001, 0.2563479777846623], "hit": 0.0},
{"timestamp": "03-27-19:25:09", "target": "jak2", "validity": 1.0, "uniqueness": 1.0, "novelty": 0.391304347826087, "top_ds": [8.675, 0.7186296483088986], "hit": 0.012422360248447204},
{"timestamp": "03-27-19:25:09", "target": "jak2", "validity": 1.0, "uniqueness": 1.0, "novelty": 0.38509316770186336, "top_ds": [8.2375, 0.3292307050426149], "hit": 0.0},
{"timestamp": "03-27-19:25:09", "target": "jak2", "validity": 1.0, "uniqueness": 1.0, "novelty": 0.34810126582278483, "top_ds": [8.328571428571427, 0.24299715851758236], "hit": 0.0}
    ]
    evalkeys = ["hit","top_ds","novelty"]
    autostat(perf,evalkeys)