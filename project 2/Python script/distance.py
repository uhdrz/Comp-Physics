import numpy as np
import matplotlib.pyplot as plt
from cmath import sqrt
from scipy.stats import gaussian_kde


x=np.genfromtxt('pos001.txt',dtype=float,usecols=(0))
y=np.genfromtxt('pos001.txt',dtype=float,usecols=(1))


numberofelectrons=2
n=0
sumdis=0

while n<len(x):
    disx=x[n]- x[n+1]
    disy=y[n]- y[n+1]
    sumdis +=sqrt(disx*disx+disy*disy)
    n=n+2
    

meandistance=sumdis/(len(x)*0.5)


print meandistance
plt.hist2d(x, y, (100, 100), cmap=plt.cm.jet)
plt.colorbar()
plt.show()

