import pandas as pd
import numpy as np
import math
from sympy import FiniteSet

import argparse
import os

parser = argparse.ArgumentParser(description='Scan statistic to identify kataegis')

args = parser.parse_args()
#print(args)

input_dir = os.path.join('data')
input_filenames = list(filter(lambda x: x.endswith(".tsv"), os.listdir(input_dir)))

#print(input_filenames)
# df = pd.read_csv()

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

