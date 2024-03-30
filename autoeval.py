from evaluate import evaluatesmis
import json


def autoeval(smiles_list,fold=4):
    for smile_path in smiles_list:
        if "parp1" in smile_path:
            evaluatesmis("parp1",smile_path,fold)
        elif "fa7" in smile_path:
            evaluatesmis("fa7",smile_path,fold)
        elif "5ht1b" in smile_path:
            evaluatesmis("5ht1b",smile_path,fold)
        elif "braf" in smile_path:
            evaluatesmis("braf",smile_path,fold)
        elif "jak2" in smile_path:
            evaluatesmis("jak2",smile_path,fold)
        else:
            raise ValueError

def loganalysis(logpath,filter_dict,valid_thres,sort_key):
    alllogs = []
    with open(logpath,"r") as f:
        for line in f:
            line = line.strip()
            resdict = json.loads(line)
            alllogs.append(resdict)
    selected = []
    for rec in alllogs:
        flag = True
        for k,v in filter_dict.items():
            if rec[k]!=v:
                flag = False
        if rec["VALID"]<valid_thres:
            flag=False
        if flag:
            selected.append(rec)
    final_list = sorted(selected,key=lambda x:x[sort_key],reverse=True)
    for rec in final_list[:5]:
        print("round: ",rec["valround"],"valid: ", rec["VALID"],"seed:",rec["seed"], "hit:", rec["hit"],"train_method:", rec["train_method"],"target_prop: ",rec["target_prop"])
    # print(final_list[:5])

# if __name__ == "__main__":
#     logpath = r"/root/MOOD/evaluation_dictzinc.log"
#     filter_dict = {
#         "train_method": "isgdpo",
#         "target_prop": "jak2",
#     }
#     valid_thres = 60
#     sort_key = "hit"
#     loganalysis(logpath,filter_dict,valid_thres,sort_key)
    

            
            

    

if __name__ == "__main__":
    smiles_list = [
        r"/root/MOOD/gensamples/ismethod/smiles_num1024_prop5ht1b_seed14730_val56.35.txt",
        r"/root/MOOD/gensamples/ismethod/smiles_num1024_prop5ht1b_seed28673_val77.73.txt",
        r"/root/MOOD/gensamples/ismethod/smiles_num1024_propbraf_seed5988_val57.91.txt",
        r"/root/MOOD/gensamples/ismethod/smiles_num1024_propbraf_seed12984_val83.01.txt",
        r"/root/MOOD/gensamples/ismethod/smiles_num1024_propfa7_seed22096_val59.08.txt",
        r"/root/MOOD/gensamples/ismethod/smiles_num1024_propfa7_seed25066_val82.32.txt",
        r"/root/MOOD/gensamples/ismethod/smiles_num1024_propjak2_seed6803_val83.98.txt",
        r"/root/MOOD/gensamples/ismethod/smiles_num1024_propjak2_seed20406_val62.40.txt",
        r"/root/MOOD/gensamples/ismethod/smiles_num1024_propparp1_seed21_val82.71.txt",
        r"/root/MOOD/gensamples/ismethod/smiles_num1024_propparp1_seed121_val91.99.txt"
    ]
    autoeval(smiles_list,4)
