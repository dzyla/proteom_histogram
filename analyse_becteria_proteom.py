from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
import pandas as pd

data_file = []
value = 15
value1 = 19

for i,rec in enumerate(SeqIO.parse('Ecoli_k12_proteom.fasta', 'fasta')):
    if 'X' not in rec.seq:
        protein = ProteinAnalysis(str(rec.seq))
        sum = 0
        distribution = protein.get_amino_acids_percent()
        for aa in ['A','I','L','M','F','V','P','W']:
            if aa in distribution:
                sum += distribution[aa]
        record = [i, round(protein.molecular_weight() / 1000, 4), round(protein.isoelectric_point(), 2),
                  sum * 100]
        data_file.append(np.array(record))
        '''optional if you want to screen for features of the protein'''
        # if name in rec.description:
        #     print(rec.description,record,rec.seq)
        # if value < record[1] < value1:
        #     print(rec.description, record)

data_file = np.asarray(data_file)
df = pd.DataFrame.from_records(data_file)
#
plt.figure(figsize=(14,8))
fig = plt.gcf()
fig.suptitle("E. coli K12 protein analysis", fontsize=14)



plt.subplot(131)
plt.hist(df[1],100,histtype='bar')
plt.xlabel('Mass distribution of E. coli K12 proteins, KDa')
plt.ylabel('Number of proteins')
#plt.xscale('log',basex=10)


plt.subplot(132)
plt.hist(df[2],80,histtype='bar',facecolor='g')
plt.xlabel('pI distribution of E. coli K12 proteins')
plt.ylabel('Number of proteins')

plt.subplot(133)
plt.hist(df[3],50,histtype='bar',facecolor='orchid')
plt.xlabel('Hydrophobic residues content [ALIMFVPW], %')
plt.ylabel('Number of proteins')

plt.show()