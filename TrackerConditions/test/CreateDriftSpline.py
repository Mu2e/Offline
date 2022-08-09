import GeneralUtilities
import numpy as np
import sys

data = np.genfromtxt(sys.argv[1])
fout = open(sys.argv[2],"w")

is_res = int(sys.argv[3]) == 1

x = data[:,0]
y = data[:,1]

if not is_res:
  # FIXME correction for negative velocities
  for j in range(len(y)-1):
    if y[j] > y[j+1]:
      y[j] = y[j+1] - 1e-5
if is_res:
  s = GeneralUtilities.SplineInterpolation(x,y,False,True)
else:
  s = GeneralUtilities.SplineInterpolation(x,y,True,True)
splineA = list(s.getSplineA())
splineB = list(s.getSplineB())
splineC = list(s.getSplineC())
splineD = list(s.getSplineD())

name = "driftSpline"
if is_res:
  name = "driftResSpline"
fout.write("services.ProditionsService.strawResponse.%sBins : [" % name)
for i in range(len(x)):
  fout.write("%.8f" % x[i])
  if i != len(x)-1:
    fout.write(", ")
fout.write("]\n")

fout.write("services.ProditionsService.strawResponse.%sA : [" % name)
for i in range(len(splineA)):
  fout.write("%.8f" % splineA[i])
  if i != len(splineA)-1:
    fout.write(", ")
fout.write("]\n")

fout.write("services.ProditionsService.strawResponse.%sB : [" % name)
for i in range(len(splineB)):
  fout.write("%.8f" % splineB[i])
  if i != len(splineB)-1:
    fout.write(", ")
fout.write("]\n")

fout.write("services.ProditionsService.strawResponse.%sC : [" % name)
for i in range(len(splineC)):
  fout.write("%.8f" % splineC[i])
  if i != len(splineC)-1:
    fout.write(", ")
fout.write("]\n")

fout.write("services.ProditionsService.strawResponse.%sD : [" % name)
for i in range(len(splineD)):
  fout.write("%.8f" % splineD[i])
  if i != len(splineD)-1:
    fout.write(", ")
fout.write("]\n")

fout.close()
