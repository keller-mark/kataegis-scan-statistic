import unittest
import os
import pandas as pd

from scan_stat import BernoulliScanStatistic, PoissonScanStatistic

class Test(unittest.TestCase):

  @classmethod
  def setUpClass(cls):

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

    cls.df = load_df()

    cls.pss = PoissonScanStatistic(cls.df)
    cls.bss = BernoulliScanStatistic(cls.df)

  def test_init(self):
    self.assertEqual(Test.pss.n_G, 249250621)
    self.assertEqual(Test.bss.n_G, 249250621)

    self.assertAlmostEqual(Test.pss.mu_G, float(4.278283683252548e-05))
    self.assertAlmostEqual(Test.bss.mu_G, float(4.278283683252548e-05))
  
  def test_monte_carlo_bernoulli(self):
    Test.bss.monte_carlo()
    Test.pss.monte_carlo()


if __name__ == '__main__':
  unittest.main()
