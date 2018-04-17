import pandas as pd

import argparse
import os
from scan_stat import BernoulliScanStatistic, PoissonScanStatistic

parser = argparse.ArgumentParser(description='Scan statistic to identify kataegis')

args = parser.parse_args()

input_dir = 'data'
input_filenames = list(filter(lambda x: x.endswith(".tsv"), os.listdir(input_dir)))
df = None
for input_filename in input_filenames:
    file_df = pd.read_csv(os.path.join(input_dir, input_filename), sep='\t')
    if df is None:
        df = file_df
    else:
        df = df.append(file_df)

bss = BernoulliScanStatistic(df)
bss.monte_carlo()

pss = PoissonScanStatistic(df)
pss.monte_carlo()

