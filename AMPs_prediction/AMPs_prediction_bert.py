from os import environ
from sys import argv

# usage python prediction_bert.py EUK_PEP_DB_drp.fa EUK.proba.tsv
environ["CUDA_VISIBLE_DEVICES"] = "1"  # 指定所要使用的显卡
seq_path = 'ASM584v2.orf.fa'
from bert_sklearn import BertClassifier
from bert_sklearn import load_model
import numpy as np
import pandas as pd

model = load_model("bert.bin")
tmp = pd.read_csv(seq_path, sep="\t", header=None, names=["seq"], index_col=False).seq.values
seq_array = []
for eachseq in tmp:
    if ">" not in eachseq:
        seq_array.append(" ".join(list(eachseq)))

seq_array = np.array(seq_array)
y_prob = model.predict_proba(seq_array)
y_prob = y_prob[:, 1]
pd.DataFrame(y_prob).to_csv('bert.txt', sep="\t", header=False, index=False)