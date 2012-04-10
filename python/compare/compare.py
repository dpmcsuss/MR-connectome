
import argparse
import matplotlib.pyplot
import numpy as np
from scipy.io import loadmat
from scipy.sparse import csc_matrix
import csv


parser = argparse.ArgumentParser(description='False color graph from a mat file')
parser.add_argument('matfile', action="store")
parser.add_argument('csvfile', action="store")

result = parser.parse_args()

matcontents = loadmat ( result.matfile )

matcsc = matcontents["fibergraph"] 
matdata = np.array (matcsc.todense())

reader = csv.reader(open(result.csvfile,"rb"))
csvdata = np.array([[float(col) for col in row] for row in reader])

for j in range (matdata.shape[1]):
  for i in range (matdata.shape[0]):
    if matdata [ i,j ] != csvdata [i+1,j+1]:
      print "Data differ at (%d,%d): values mat %d csv %d" % (i,j,matdata[i,j],csvdata[i,j])



