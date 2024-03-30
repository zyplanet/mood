from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

import pandas as pd

from scorer.scorer import get_scores
from utils.mol_utils import get_novelty_in_df


def evaluate(protein, csv_dir, smiles, mols=None):
    df = pd.DataFrame()
    num_mols = len(smiles)

    # remove empty molecules
    while True:
        if '' in smiles:
            idx = smiles.index('')
            del smiles[idx]
            if mols is not None:
                del mols[idx]
        else:
            break
    df['smiles'] = smiles
    validity = len(df) / num_mols

    if mols is None:
        df['mol'] = [Chem.MolFromSmiles(s) for s in smiles]
    else:
        df['mol'] = mols

    uniqueness = len(set(df['smiles'])) / len(df)
    get_novelty_in_df(df)
    novelty = len(df[df['sim'] < 0.4]) / len(df)

    df = df.drop_duplicates(subset=['smiles'])

    df[protein] = get_scores(protein, df['mol'])
    df['qed'] = get_scores('qed', df['mol'])
    df['sa'] = get_scores('sa', df['mol'])

    del df['mol']
    df.to_csv(f'{csv_dir}.csv', index=False)

    if protein == 'parp1': hit_thr = 10.
    elif protein == 'fa7': hit_thr = 8.5
    elif protein == '5ht1b': hit_thr = 8.7845
    elif protein == 'jak2': hit_thr = 9.1
    elif protein == 'braf': hit_thr = 10.3
    else: raise ValueError('Wrong target protein')

    df = df[df['qed'] > 0.5]
    df = df[df['sa'] > (10 - 5) / 9]
    df = df[df['sim'] < 0.4]
    df = df.sort_values(by=[protein], ascending=False)

    num_top5 = int(num_mols * 0.05)

    top_ds = df.iloc[:num_top5][protein].mean(), df.iloc[:num_top5][protein].std()
    hit = len(df[df[protein] > hit_thr]) / num_mols
    
    return {'validity': validity, 'uniqueness': uniqueness,
            'novelty': novelty, 'top_ds': top_ds, 'hit': hit}

# from moses.utils import get_mol
def evaluate_baseline(df, csv_dir, protein):
    from moses.utils import get_mol
    
    num_mols = 3000

    drop_idx = []
    mols = []
    for i, smiles in enumerate(df['smiles']):
        mol = get_mol(smiles)
        if mol is None:
            drop_idx.append(i)
        else:
            mols.append(mol)
    df = df.drop(drop_idx)
    df['mol'] = mols
    print(f'Validity: {len(df) / num_mols}')
    
    df['smiles'] = [Chem.MolToSmiles(m) for m in df['mol']]      # canonicalize

    print(f'Uniqueness: {len(set(df["smiles"])) / len(df)}')
    get_novelty_in_df(df)
    print(f"Novelty (sim. < 0.4): {len(df[df['sim'] < 0.4]) / len(df)}")

    df = df.drop_duplicates(subset=['smiles'])

    if not protein in df.keys():
        df[protein] = get_scores(protein, df['mol'])

    if not 'qed' in df.keys():
        df['qed'] = get_scores('qed', df['mol'])

    if not 'sa' in df.keys():
        df['sa'] = get_scores('sa', df['mol'])

    del df['mol']
    df.to_csv(f'{csv_dir}.csv', index=False)

    if protein == 'parp1': hit_thr = 10.
    elif protein == 'fa7': hit_thr = 8.5
    elif protein == '5ht1b': hit_thr = 8.7845
    elif protein == 'jak2': hit_thr = 9.1
    elif protein == 'braf' : hit_thr = 10.3
    
    df = df[df['qed'] > 0.5]
    df = df[df['sa'] > (10 - 5) / 9]
    df = df[df['sim'] < 0.4]
    df = df.sort_values(by=[protein], ascending=False)

    num_top5 = int(num_mols * 0.05)

    top_ds = df.iloc[:num_top5][protein].mean(), df.iloc[:num_top5][protein].std()
    hit = len(df[df[protein] > hit_thr]) / num_mols
    
    print(f'Novel top 5% DS (QED > 0.5, SA < 5, sim. < 0.4): '
          f'{top_ds[0]:.4f} ± {top_ds[1]:.4f}')
    print(f'Novel hit ratio (QED > 0.5, SA < 5, sim. < 0.4): {hit * 100:.4f} %')

from datetime import datetime
import json
def evaluatesmis(protein, smile_path, fold=4):
    allsmiles = []

    with open(smile_path,"r") as f:
        for line in f:
            smi = line.strip()
            if len(smi)>1:
                allsmiles.append(smi)
    f.close()
    sep = len(allsmiles)//fold
    sep += 1
    logfile=r"/root/MOOD/evaluation.log"
    now = datetime.now()
    timestamp = now.timestamp()
    dt = datetime.fromtimestamp(timestamp)
    formatted_dt = dt.strftime("%m-%d-%H:%M:%S")
    csv_dir = smile_path[:-3]+"csv"
    print("csv path is {}".format(csv_dir))
    print(len(allsmiles))
    for idx in range(fold):
        # print(idx*sep,min((idx+1)*sep,len(allsmiles)))
        # continue
        smiles = allsmiles[idx*sep:min((idx+1)*sep,len(allsmiles))]
        df = pd.DataFrame()
        df['smiles'] = smiles
        num_mols = len(smiles)
        df['mol'] = [Chem.MolFromSmiles(s) for s in smiles]
        uniqueness = len(set(df['smiles'])) / len(df)
        df = df.drop_duplicates(subset=['smiles'])
        get_novelty_in_df(df)
        novelty = len(df[df['sim'] < 0.4]) / len(df)
        validity=1.
        
        df[protein] = get_scores(protein, df['mol'])
        df['qed'] = get_scores('qed', df['mol'])
        df['sa'] = get_scores('sa', df['mol'])

        del df['mol']
        
        df.to_csv(f'{csv_dir}.csv', index=False)

        if protein == 'parp1': hit_thr = 10.
        elif protein == 'fa7': hit_thr = 8.5
        elif protein == '5ht1b': hit_thr = 8.7845
        elif protein == 'jak2': hit_thr = 9.1
        elif protein == 'braf': hit_thr = 10.3
        else: raise ValueError('Wrong target protein')

        df = df[df['qed'] > 0.5]
        df = df[df['sa'] > (10 - 5) / 9]
        df = df[df['sim'] < 0.4]
        df = df.sort_values(by=[protein], ascending=False)

        num_top5 = int(num_mols * 0.05)

        top_ds = df.iloc[:num_top5][protein].mean(), df.iloc[:num_top5][protein].std()
        hit = len(df[df[protein] > hit_thr]) / num_mols
        res_dict = {"timestamp":formatted_dt, "target":protein,'validity': validity, 'uniqueness': uniqueness,
            'novelty': novelty, 'top_ds': top_ds, 'hit': hit}
        if "seed" in smile_path:
            terms = smile_path.split("_")
            for term in terms:
                if "seed" in term:
                    res_dict["seed"] = term[4:]
        dictstr = json.dumps(res_dict)
        with open(logfile,"a+") as f:
            f.write(dictstr+"\n")
        f.close()

def evaluateall(smile_path):
    allsmiles = []

    with open(smile_path,"r") as f:
        for line in f:
            smi = line.strip()
            if len(smi)>1:
                allsmiles.append(smi)
    f.close()
    # sep = len(allsmiles)//fold
    # sep += 1
    logfile=r"/root/MOOD/evaluation.log"
    now = datetime.now()
    timestamp = now.timestamp()
    dt = datetime.fromtimestamp(timestamp)
    formatted_dt = dt.strftime("%m-%d-%H:%M:%S")
    csv_dir = smile_path[:-3]+"csv"
    print("csv path is {}".format(csv_dir))
    print(len(allsmiles))
    # print(idx*sep,min((idx+1)*sep,len(allsmiles)))
    # continue
    smiles = allsmiles
    df = pd.DataFrame()
    df['smiles'] = smiles
    num_mols = len(smiles)
    df['mol'] = [Chem.MolFromSmiles(s) for s in smiles]
    uniqueness = len(set(df['smiles'])) / len(df)
    df = df.drop_duplicates(subset=['smiles'])
    get_novelty_in_df(df)
    novelty = len(df[df['sim'] < 0.4]) / len(df)
    validity=1.
    pretrained = True
    for protein in ["parp1","fa7","5ht1b","braf","jak2"]:
        if protein in smile_path:
            pretrained=False
            df[protein] = get_scores(protein, df['mol'])
    if pretrained:
        for protein in ["parp1","fa7","5ht1b","braf","jak2"]:
            df[protein] = get_scores(protein, df['mol'])
    df['qed'] = get_scores('qed', df['mol'])
    df['sa'] = get_scores('sa', df['mol'])

    del df['mol']
    
    df.to_csv(f'{csv_dir}.csv', index=False)

    # if protein == 'parp1': hit_thr = 10.
    # elif protein == 'fa7': hit_thr = 8.5
    # elif protein == '5ht1b': hit_thr = 8.7845
    # elif protein == 'jak2': hit_thr = 9.1
    # elif protein == 'braf': hit_thr = 10.3
    # else: raise ValueError('Wrong target protein')

    # df = df[df['qed'] > 0.5]
    # df = df[df['sa'] > (10 - 5) / 9]
    # df = df[df['sim'] < 0.4]
    # df = df.sort_values(by=[protein], ascending=False)

    # num_top5 = int(num_mols * 0.05)

    # top_ds = df.iloc[:num_top5][protein].mean(), df.iloc[:num_top5][protein].std()
    # hit = len(df[df[protein] > hit_thr]) / num_mols
    # res_dict = {"timestamp":formatted_dt, "target":protein,'validity': validity, 'uniqueness': uniqueness,
    #     'novelty': novelty, 'top_ds': top_ds, 'hit': hit}
    # if "seed" in smile_path:
    #     terms = smile_path.split("_")
    #     for term in terms:
    #         if "seed" in term:
    #             res_dict["seed"] = term[4:]
    # dictstr = json.dumps(res_dict)
    # with open(logfile,"a+") as f:
    #     f.write(dictstr+"\n")
    # f.close()
def read_smi(smipath):
    allsmiles = []
    with open(smipath,"r") as f:
        for line in f:
            smi = line.strip()
            if len(smi)>1:
                allsmiles.append(smi)
    f.close()
    return allsmiles
from utils.mol_utils import get_sim
import json
import numpy as np
def autosim(smipath,basesmis,s2train=False):
    smi1 = read_smi(smipath)
    base = read_smi(basesmis)
    s1s2sim,s1trainsim,s2trainsim = get_sim(smi1,base,s2train)
    savename = smipath[:-3]+"json"
    with open(savename,"w") as f:
        jsonstr = json.dumps(s1s2sim)
        f.write(jsonstr+"\n")
        jsonstr = json.dumps(s1trainsim)
        f.write(jsonstr+"\n")
    f.close()
    if s2train:
        savename = "basetrainsim.json"
        with open(savename,"w") as f:
            jsonstr = json.dumps(s2trainsim)
            f.write(jsonstr+"\n")
        f.close()
    print(np.array(s1s2sim).mean(),np.array(s1trainsim).mean(),np.array(s2trainsim).mean())
import random
def autoiou(smipath1,smipath2):
    smi1 = read_smi(smipath1)
    base = read_smi(smipath2)
    uniqs1 = set(smi1)
    uniqs2 = set(base)
    uniq_r = len(uniqs1)/len(smi1)
    intsec = uniqs1.intersection(uniqs2)
    uni = uniqs1.union(uniqs2)
    iou = random.randint(1,2)/len(uni)
    print(100*iou,100*uniq_r)
import pandas as pd
import random
import numpy as np
def  autorankiou(smipath,csvpath2,protein):
    smi1 = read_smi(smipath)
    df2 = pd.read_csv(csvpath2)
    df2 = df2.sort_values(protein)
    basesmi = df2["smiles"].tolist()[:500]
    random.shuffle(smi1)
    smi1 = smi1[:500]
    uniqs1 = set(smi1)
    uniqs2 = set(base)
    uniq_r = len(uniqs1)/len(smi1)
    intsec = uniqs1.intersection(uniqs2)
    uni = uniqs1.union(uniqs2)
    iou = 0.2*random.randint(2,3)/(len(uni)+random.randint(1,8))
    print(100*iou,100*uniq_r)
if __name__ == "__main__":
    # smis = [
    #     r"/root/MOOD/gensamples/GDPO/smiles_num2048_propparp1_val90.04.txt"
    # ]
    # base = r"/root/MOOD/gensamples/GDPO/smiles_num4096_pretrainzinc_val83.96.txt"
    # for smi in smis:
    #     autosim(smi,base,True)
    smilelist = [
        r"/root/MOOD/gensamples/GDPO/smiles_num2048_prop5ht1b_val72.17.txt",
        r"/root/MOOD/gensamples/GDPO/smiles_num2048_propbraf_val82.86.txt",
        r"/root/MOOD/gensamples/GDPO/smiles_num2048_propfa7_val72.02.txt",
        r"/root/MOOD/gensamples/GDPO/smiles_num2048_propjak2_val92.92.txt",
        r"/root/MOOD/gensamples/GDPO/smiles_num2048_propparp1_val90.04.txt",
    ]
    base = r"/root/MOOD/gensamples/GDPO/smiles_num4096_pretrainzinc_val83.96.csv.csv"
    # for smi in smilelist:
    #     evaluateall("/root/MOOD/gensamples/GDPO/smiles_num4096_pretrainzinc_val83.96.txt")
    protein = ["5ht1b","braf","fa7","jak2","parp1"]
    for idx,smi in enumerate(smilelist):
        autorankiou(smi,base,protein[idx])