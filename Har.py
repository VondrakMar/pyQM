import numpy as np

def STO_H(r,R1):
    return 0.3696*np.exp(-0.4166*np.linalg.norm(r-R1)**2)
def STO_He(r,R1):
    return 0.5881*np.exp(-0.7739*np.linalg.norm(r-R1)**2)



import matplotlib.pyplot as plt


r1 = np.array([0.0,0.0,0.0])
R1s = np.linspace(0,3,1000)
STO_Hvalues = []
STO_Hevalues = []
for R1 in R1s:
    STO_Hvalues.append(STO_H(r1,R1))
    STO_Hevalues.append(STO_He(r1,R1))
plt.plot(R1s,STO_Hvalues)
plt.plot(R1s,STO_Hevalues)
plt.show()

### integrate potential for hydrogen atom
def Vrs(phi1,phi2,Z,r):
    
