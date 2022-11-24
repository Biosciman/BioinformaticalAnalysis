from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

P1 = 'effector_exp_protein.fasta'
P2 = 'T6SE_Training_Pos_138.fasta'
N1 = 'T6SE_Training_Neg_1112.fasta'

Gene_original = []
Gene_seq = []
for record in SeqIO.parse(P1, "fasta"):
    b = str(record.seq.strip('*'))
    if "*" not in b:
        Gene_original.append([record.id, record.description, b])
        Gene_seq.append(b)

a = 0
for record in SeqIO.parse(P2, "fasta"):
    b = str(record.seq.strip())
    if "*" not in b:
        if b in Gene_seq:
            a += 1
            # print('The protein is exist! Number:', a)
        else:
            # print('OK!')
            Gene_seq.append(str(record.seq.strip()))
            Gene_original.append([record.id, record.description, str(record.seq.strip())])

# 返回list中最长字符串
# res = max(Gene_seq, key=len, default='') # 1623个氨基酸
# print(len(res))

# with open('T6SS_Negative_length.csv', 'w') as f:
#     for record in SeqIO.parse(N1, "fasta"):
#         f.write(str(len(str(record.seq.strip('')))))
#         f.write('\n')


with open('T6SS_Positive_length.csv', 'w') as f:
    for i in Gene_seq:
        f.write(str(len(i)))
        f.write('\n')


# 写入新的fasta文件
rec = []
for i in Gene_original:
        rec.append(SeqRecord(
            Seq(i[2]),
            id=i[0],
            description=i[1]
        ))

SeqIO.write(rec, 'T6SS_Positive.fa', 'fasta')




# Gene_original = []
# Gene_seq = []
# for record in SeqIO.parse(N1, "fasta"):
#     Gene_original.append([record.id, record.description, str(record.seq.strip(''))])
#     Gene_seq.append(str(record.seq.strip('')))
#
# # 写入新的fasta文件
# rec = []
# for i in Gene_original:
#     # print(len(str(record.seq.strip(''))))
#     if len(str(record.seq.strip(''))) < 2000:
#         rec.append(SeqRecord(
#             Seq(i[2]),
#             id=i[0],
#             description=i[1]
#         ))
# SeqIO.write(rec, 'T6SS_Negative.fa', 'fasta')