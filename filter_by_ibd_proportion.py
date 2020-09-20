#! /usr/bin/env python3

import pandas as pd
import sys

if len(sys.argv) != 4:
    print("usage:\n"
          "  ./filter_by_ibd_proportion.py xxxx.ibd <int_genome_size_cM> <dbl_ratio>\n\n")
    sys.exit(2)
ibd_file = sys.argv[1]
genome_size = int(sys.argv[2])
threshold_ratio = float(sys.argv[3])

# read in 
df = pd.read_csv(ibd_file, header=None, sep='\t')
df.columns = ['id1', 'hap1', 'id2', 'hap2', 'chr', 'start', 'end', 'score', 'cM']

# calculate ibdtotal
df2 = df.groupby(['id1', 'id2'])['cM'].sum()
df2 = df2.reset_index()
df2 = df2.rename(columns={'cM': 'ibdtotal'})

# add ibdtotal to df
df = df.merge(df2, on = ['id1', 'id2'], how='inner')

# filter by ibdtotal
df = df[df['ibdtotal'] < genome_size * threshold_ratio]
print(df.shape, file=sys.stderr)

# write out to stdout
df.iloc[:, 0:-1].to_csv(sys.stdout, index=False, header=None, sep='\t')
