import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description='Scan statistic to identify kataegis')

args = parser.parse_args()
print(args)

input_dir = os.path.join('data')
input_filenames = list(filter(lambda x: x.endswith(".tsv"), os.listdir(input_dir)))

print(input_filenames)
# df = pd.read_csv()
