from collections import defaultdict

# 此为取消0 padding的版本。

seq_name_dict = defaultdict(str)
with open('T6SS_Negative.fa') as f:
	for line in f:
		line = line.rstrip()
		if line.startswith(">"):
			name = line[1:]
		else:
			seq_name_dict[name] += line
# print(seq_name_dict)
alphabet = "ACDEFGHIKLMNPQRSTVWY"
aa_dict = {}
for i, aa in enumerate(alphabet):
	aa_dict[aa] = i+1
#print(aa, aa_dict[aa], '\n')
#

sequence = []
for i in seq_name_dict.values():
	a = ""
	count = 0
	for letter in i:
		b = aa_dict[letter]
		a += str(b)
		count += 1
		# print(count)
		# print(len(i))
		if count != len(i):
			a += ','
		else:
			sequence.append(a)
print(sequence)

with open('T6SS_Negative2.txt','w') as f:
	for i in sequence:
		f.write(i)
		f.write('\n')
