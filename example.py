import os
import pandas as pd

from scan_stat import BernoulliScanStatistic, PoissonScanStatistic

def load_df(input_dir = 'data'):
    input_filenames = list(filter(lambda x: x.endswith(".tsv"), os.listdir(input_dir)))
    df = None
    for input_filename in input_filenames:
        file_df = pd.read_csv(os.path.join(input_dir, input_filename), sep='\t')
        if df is None:
            df = file_df
        else:
            df = df.append(file_df)
    return df

if __name__ == '__main__':
    df = load_df()
    bss = BernoulliScanStatistic(df)
    bss.run()