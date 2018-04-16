import numpy as np
from scipy.special import comb
import math
from flanking import Flanking


class BernoulliScanStatistic():

  """
    Reference: https://www.satscan.org/papers/k-cstm1997.pdf

    H_0 : p = q
    H_1 : p > q, zone Z in window

    mu(A) = integer for all subsets A of space G


    n_Z = observed number of points in zone Z
    n_G = total number of observed points

    point = individual in a state of kataegis or not

    p = probability that each individual within the zone is a point
    q = probability that each individual outside the zone is a point

  """

  def __init__(self, mutation_df):
    self.flanking = Flanking('hg')
    self.chromosomes = self.flanking.get_chromosome_lengths()
    print(self.chromosomes)
    self.chromosome = self.chromosomes.keys()[0]
    print(self.chromosome)
    
    self.mutation_df = mutation_df
    self.n_G = self.chromosomes[self.chromosome]
    self.mu_G = len(mutation_df.loc["chr" + mutation_df['Chromosome'] == self.chromosome])
    print(self.mu_G)


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

  
  """
  Find the distribution of the test statistic using Monte Carlo simulation

  Using underlying measure mu, obtain replications of the data set generated under H_0,
    conditioning on total number of points n_G.

  With 9999 replications, test is significant at alpha = 0.05
    if test statistic for real data set is among 500 highest values 
    of the test statistic from the replications

  """
  def monte_carlo(self):
    n_G = self.n_G
    mu_G = self.mu_G

    def binomial_pmf(n, p, k):
      # n = number of trials
      # p = success probability in each trial
      # k = number of successes
      return comb(n, k) * np.power(p, k) * np.power((1 - p), (n - k))

    n_trials = 9999
    alpha = 0.05
    n_top = math.ceil(n_trials * alpha)

    for _ in range(0, n_trials):
      print(binomial_pmf(n_G, mu_G, ))
    

  




