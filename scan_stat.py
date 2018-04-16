import numpy as np
import math
from sympy import FiniteSet


class BernoulliScanStatistic():

  """

    H_0 : p = q
    H_1 : p > q, zone Z in window

    mu(A) = integer for all subsets A of space G


    n_Z = observed number of points in zone Z
    n_G = total number of observed points

    point = individual in a state of kataegis or not

    p = probability that each individual within the zone is a point
    q = probability that each individual outside the zone is a point

  """

  def __init__(self):
    pass


  """
  Probability of n_G number of points in the study area
  """
  def prob_n_G(n_G, p, q, mu_G, mu_Z):
    return (
        np.exp( -p*mu_Z - q*(mu_G - mu_Z) )
      * np.power(
          (p*mu_Z + q*(mu_G - mu_Z)), 
          n_G
        )
    ) / math.factorial(n_G)

  """
  pdf f(x) of specific point being observed at location x
  """
  def pdf(x, p, q, mu_G, mu_Z, x_in_Z):
    mu_x = 0 # TODO: ???
    if x_in_Z:
      return (p*mu_x) / (p*mu_Z + q*(mu_G - mu_Z))
    else:
      return (q*mu_x) / (p*mu_Z + q*(mu_G - mu_Z))

  """
  Likelihood function L(Z, p, q)
  """
  def likelihood(Z, p, q, n_G, mu_G, mu_Z):
    xs_in_Z = [] # TODO: array of x's in Z
    xs_not_in_Z = [] # TODO: array of x's not in Z

    return (
        prob_n_G(n_G, p, q, mu_G, mu_Z)
      * np.prod(
          np.apply_along_axis(
            pdf, # function
            0, # axis
            xs_in_Z, # array
            p, q, mu_G, mu_Z, True # additional args
          ) 
        )
      * np.prod(
          np.apply_along_axis(
            pdf, # function
            0, # axis
            xs_not_in_Z, # array
            p, q, mu_G, mu_Z, False # additional args
          ) 
        )
    )

  def test_statistic(Z):
    supremum = FiniteSet(*Z).sup
    return (
      (supremum * likelihood(Z, p, q, n_G, mu_G, mu_Z)) /
      (
        (np.exp(-n_G) / math.factorial(n_G)) 
      * np.power(n_G / mu_G, n_G)
      * np.prod(
          np.apply_along_axis(
            np.mean,
            0,
            Z
          ) 
        )
      )
    )

