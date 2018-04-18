import numpy as np
import math
from hg import HG

POS = 'Chromosome Start'
CHR = 'Chromosome'
DONOR_ID = 'Donor ID'

class ScanStatistic():

  """
    Motivation:
    "...Since kataegis is observed at different rates in different cancer types, 
    a tool that can dynamically optimize kataegis detection per cancer type
    and assess significance of kataegis events would be of a great use..."

    "...uses a sliding window (of fixed width) approach to test deviation of 
    observed SNV trinucleotide content and inter-mutational distance from expected 
    by chance alone." 


    Scan Statistic Reference: 
    https://www.satscan.org/papers/k-cstm1997.pdf

    H_0 : p = q
    H_1 : p > q, zone Z in window

    mu(A) = mean for all subsets A of space G

    n_Z = observed number of points in zone Z
    n_G = total number of observed points

    point = mutated base pair or not

    p = probability that each individual within the zone is a point
    q = probability that each individual outside the zone is a point

  """

  def __init__(self, mut_df):
    hg = HG('hg')
    
    chromosome = "1"
    self.chr_len = hg.get_chromosome_length("chr" + chromosome)

    self.mut_df = mut_df.loc[mut_df[CHR] == chromosome]
    print(self.mut_df)

    # total number of base pairs for chromosome
    self.n_G = self.chr_len

    self.n_mutations = len(mut_df)
    self.n_donors = len(mut_df[DONOR_ID].unique())
    # mean mutation frequency for chromosome
    self.mu_G = (self.n_mutations / self.n_donors) / self.n_G

    self.L_0 = self.likelihood_null()
  
  def get_window(self, size, offset, donor_id='DO51965'):
    min_pos = offset
    max_pos = offset + size
    window_df = self.mut_df.loc[(self.mut_df[POS] >= min_pos) & (self.mut_df[POS] < max_pos)]
    n_window_mutations = len(window_df)

    # TODO: use np.convolve to compute these means at the same time for all windows
    p = (n_window_mutations / self.n_donors) / size
    q = ((self.n_mutations - n_window_mutations) / self.n_donors) / (self.chr_len - size)

    n_window_mutations_donor = len(window_df.loc[window_df[DONOR_ID] == donor_id])
    return {
      'p': p,
      'q': q,
      'n_Z': size,
      'mu_Z': (n_window_mutations_donor / size)
    }
  
  """
  Find the distribution of the test statistic
    of the likelihood ratio (lambda) using Monte Carlo simulation

  Using underlying measure mu, obtain replications of the data set generated under H_0,
    conditioning on total number of points n_G.

  With 9999 replications, test is significant at alpha = 0.05
    if test statistic for real data set is among 500 highest values 
    of the test statistic from the replications

  """
  def monte_carlo(self):
    n_G = self.n_G
    mu_G = self.mu_G

    n_trials = 9999
    alpha = 0.05
    n_top = math.ceil(n_trials * alpha)

    poisson_lambda = self.mu_G * self.n_G
    simulations = np.random.poisson(poisson_lambda, (n_trials))
    top = np.partition(simulations, n_trials - n_top)
    top_range = [np.amin(top), np.amax(top)]
    print(top_range)

class BernoulliScanStatistic(ScanStatistic):

  def __init__(self, mut_df):
    super().__init__(mut_df)

  """
  Likelihood function L(Z, p, q)
  """
  def likelihood(self, Z):
    n_Z = Z['n_Z']
    n_G = self.n_G

    mu_Z = Z['mu_Z']
    mu_G = self.mu_G

    p = Z['p']
    q = Z['q']

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
  
  def likelihood_ratio(self, Z):
    return (self.likelihood(Z) / self.L_0)
  
  def run(self):
    start = 0
    end = self.chr_len
    window_size = 1000
    for i in range(0, end, window_size):
      Z = self.get_window(size=window_size, offset=i)
      if Z['mu_Z'] > 0:
        l = self.likelihood_ratio(Z)
        print(Z, l)

class PoissonScanStatistic(ScanStatistic):

  def __init__(self, mut_df):
    super().__init__(mut_df)