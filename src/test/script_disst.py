from pylab import *

dat=loadtxt("out.txt")
hist(dat, bins=50)
show()
