import numpy as np
import matplotlib.pyplot as plt
from cmath import sqrt
from scipy.stats import gaussian_kde


x=np.genfromtxt('pos.dat',dtype=float,usecols=(0))
y=np.genfromtxt('pos.dat',dtype=float,usecols=(1))


numberofelectrons=2
n=0
sumdis=0

while n<len(x):
    disx=x[n]- x[n+1]
    disy=y[n]- y[n+1]
    sumdis +=sqrt(disx*disx+disy*disy)
    n=n+2
    

meandistance=sumdis/(len(x)*0.5)

#N=50
#colors = np.random.rand(N)
#area = np.pi * (15 * np.random.rand(N))**2

#xy = np.vstack([x,y])
#z = gaussian_kde(xy)(xy)


#plt.scatter(x, y, c=z, s=100, edgecolor='')
#plt.show()

plt.hist2d(x, y, (100, 100), cmap=plt.cm.jet)
plt.colorbar()
plt.show()
print meandistance
