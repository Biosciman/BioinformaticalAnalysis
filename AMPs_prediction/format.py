import numpy as np
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='VAE')
parser.add_argument('--path', type=str, help='Path of sequence file')
# parser.add_argument('--flag', type=str, help='Whether the sequence is AMP, non AMP, or test')
parser.add_argument('--output', type=str, help='Output text file name')
args = parser.parse_args()

seq_path = args.path
seq_name_dict = defaultdict(str)
with open(seq_path) as f:
	for line in f:
		line = line.rstrip()
		if line.startswith(">"):
			name = line[1:]
		else:
			seq_name_dict[name] += line
pad_dict = defaultdict(str)
padlen = 300
for key in seq_name_dict.keys():
	sl = len(seq_name_dict[key])
	if sl < padlen:
		pad_dict[key] = 'X' * (padlen - sl) + seq_name_dict[key]
	if sl > padlen:
		pad_dict[key] = seq_name_dict[key][0:padlen-1]
alphabet = "ACDEFGHIKLMNPQRSTVWY"
aa_dict = {}
for i, aa in enumerate(alphabet):
	aa_dict[aa] = i+1
#print(aa, aa_dict[aa], '\n')
encoding = np.zeros((len(pad_dict.keys()), padlen), dtype=int)
for i, key in enumerate(pad_dict.keys()):
	sequence = pad_dict[key]
	for j, letter in enumerate(sequence):
		if letter in aa_dict:
			encoding[i, j] = aa_dict[letter]
#print(encoding[0:2,])
# if args.flag != 'none':
# 	if args.flag == 'P':
# #label = np.ones((len(pad_dict.keys()), 1), dtype = int)
#         encoding = np.array(np.c_[encoding, np.ones(encoding.shape[0])],dtype=int)
# 	if args.flag == 'N':
# #label = np.zeros((len(pad_dict.keys()), 1), dtype = int)
#         encoding = np.array(np.c_[encoding, np.zeros(encoding.shape[0])],dtype=int)
# #encoding = np.insert(encoding, padlen, label, axis = 1)
#print(encoding[0:2,])
np.savetxt(args.output, encoding, delimiter=',')