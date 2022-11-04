#!/usr/bin/python3

from numpy import savetxt, where
from pandas import read_csv

atch1 = read_csv("Cohort1.csv", index_col = 0)
atch2 = read_csv("Cohort2.csv", index_col = 0)
atch3 = read_csv("Cohort3.csv", index_col = 0)
atch4 = read_csv("Cohort4.csv", index_col = 0)
atch5 = read_csv("Cohort5.csv", index_col = 0)
atch6 = read_csv("Cohort6.csv", index_col = 0)
atch7 = read_csv("Cohort7.csv", index_col = 0)
atch8 = read_csv("Cohort8.csv", index_col = 0)
atch9 = read_csv("Cohort9.csv", index_col = 0)
atch10 = read_csv("Cohort10.csv", index_col = 0)
atch11 = read_csv("Cohort11.csv", index_col = 0)
atch12 = read_csv("Cohort12.csv", index_col = 0)
atch13 = read_csv("Cohort13.csv", index_col = 0)
atch14 = read_csv("Cohort14.csv", index_col = 0)
atch15 = read_csv("Cohort15.csv", index_col = 0)

l_name = [atch1, atch2, atch3, atch4, atch5, atch6, atch7, atch8, atch9, atch10, atch11,  atch12, atch13, atch14, atch15]

def get_dic(at4csv):
  colname = list(at4csv.columns)
  rowname = list(at4csv.index)
  ind = where(at4csv < 0)
  row_index = ind[0].tolist()
  col_index = ind[1].tolist()
  for i in range(len(rowname)):
      rowname[i] = rowname[i].replace("-", "_")
  a_dic = {}
  for i_ind in range(len(col_index)):
      new_key = rowname[row_index[i_ind]] + "," + colname[col_index[i_ind]]
      if new_key in a_dic.keys():
          pass
      else:
          a_dic[new_key] = namestr(at4csv, globals())[0]
  return a_dic

def merge_dict(x,y):
  for k,v in x.items():
      if k in y.keys():
          y[k] += "," + v
      else:
          y[k] = v
  return y

def namestr(obj, namespace):
  return [name for name in namespace if namespace[name] is obj]

a_new = {}
for it in l_name:
  a_new = merge_dict(a_new, get_dic(it))

pr_l = ["target,soucre,numbers"]
for i_key in a_new.keys():
  num = a_new[i_key].count(",") + 1
  if num >=  int(len(l_name)/2):
      pr = i_key + "," + str(a_new[i_key])
      #pr = i_key + "," + str(num) # show number of corhorts
      pr_l.append(pr)

savetxt("Cytoscape_file_AtLeast7nameG.csv", pr_l, fmt = "%s")
