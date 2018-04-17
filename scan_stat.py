import numpy as np
import math
from hg import HG


class ScanStatistic():

  """
    Reference: https://www.satscan.org/papers/k-cstm1997.pdf

    H_0 : p = q
    H_1 : p > q, zone Z in window

    mu(A) = mean for all subsets A of space G

    n_Z = observed number of points in zone Z
    n_G = total number of observed points

    point = mutated base pair or not

    p = probability that each individual within the zone is a point
    q = probability that each individual outside the zone is a point

  """

  def __init__(self, mutation_df):
    self.hg = HG('hg')
    self.chr_lens = self.hg.get_chromosome_lengths()

    self.mutation_df = mutation_df

    def add_absolute_pos_column(df):
      chr_lens_cum = self.hg.get_chromosomes_cumulative()
      df['Absolute Position'] = df.apply(lambda row: chr_lens_cum["chr" + row['Chromosome']] + int(row['Chromosome Start']), axis=1)
      return df
    self.mutation_df = add_absolute_pos_column(self.mutation_df)

    # total number of base pairs in human genome
    self.n_G = np.sum(list(self.chr_lens.values()))

    n_mutations = len(mutation_df)
    n_donors = len(mutation_df['Donor ID'].unique())
    # mean mutation frequency
    self.mu_G = (n_mutations / n_donors) / self.n_G
  
  """
  Find the distribution of the test statistic using Monte Carlo simulation

  Using underlying measure mu, obtain replications of the data set generated under H_0,
    conditioning on total number of points n_G.

  With 9999 replications, test is significant at alpha = 0.05
    if test statistic for real data set is among 500 highest values 
    of the test statistic from the replications

  """
  def monte_carlo(self, np_rand_fun, args):
    n_G = self.n_G
    mu_G = self.mu_G

    n_trials = 9999
    alpha = 0.05
    n_top = math.ceil(n_trials * alpha)

    simulations = np_rand_fun(*args, (n_trials))
    top = np.partition(simulations, n_trials - n_top)
    top_range = [np.amin(top), np.amax(top)]
    print(top_range)

class BernoulliScanStatistic(ScanStatistic):

  def __init__(self, mutation_df):
    super().__init__(mutation_df)

  """
  Likelihood function L(Z, p, q)
  """
  def likelihood(self, Z, p, q):
    n_Z = len(Z)
    n_G = self.n_G

    mu_Z = np.mean(Z)
    mu_G = self.mu_G

    return (
      np.power(p, n_Z)
      * np.power((1 - p), mu_Z - n_Z)
      * np.power(q, n_G - n_Z)
      * np.power((1 - q), (mu_G - mu_Z) - (n_G - n_Z))
    )

  def likelihood_null(self):
    n_G = self.n_G
    mu_G = self.mu_G
    return (
      np.power((n_G / mu_G), n_G) 
      * np.power(((mu_G - n_G) / mu_G), (mu_G - n_G))
    )

  def monte_carlo(self):
    return super().monte_carlo(np.random.binomial, args=[self.n_G, self.mu_G])

class PoissonScanStatistic(ScanStatistic):

  def __init__(self, mutation_df):
    super().__init__(mutation_df)

  def monte_carlo(self):
    return super().monte_carlo(np.random.poisson, args=[self.mu_G])