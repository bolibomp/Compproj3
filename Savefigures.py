import matplotlib.pyplot as plt
import numpy as np

histoT1=[[375729, 1134377, 1482400, 1119836, 583946, 224882, 63932, 13013, 1764, 121, 0, 0], 
         [358922, 829950, 1190103, 1085936, 822470, 444640, 187888, 64445, 12023, 3364, 259, 0], 
         [1342960, 2060489, 1115576, 366164, 93856, 18125, 2554, 229, 47, 0, 0, 0], 
         [617194, 1174974, 1316902, 1016878, 568129, 227362, 63977, 12938, 1563, 78, 5, 0], 
         [2608754, 1655457, 564250, 140859, 26533, 3790, 357, 0, 0, 0, 0, 0]]
histoT10=[[1441018, 1925983, 1127511, 390504, 94206, 18037, 2454, 270, 17, 0, 0, 0], 
          [1629715, 1652047, 1064668, 445744, 156663, 40588, 8909, 1498, 151, 17, 0, 0], 
          [2657668, 1809568, 450717, 71897, 9113, 957, 71, 9, 0, 0, 0, 0], 
          [2002122, 1715766, 874520, 309674, 80115, 15339, 2227, 223, 13, 1, 0, 0], 
          [3739104, 1070136, 168478, 20347, 1828, 106, 1, 0, 0, 0, 0, 0]]
def boltzmaN(hist,T):
    G=5768299665
    u=[]
    g=[]
    for i in range(len(hist)):
        u.append(hist[i]*np.exp(-i/T))
    C=np.sum(u)    
    for i in range(len(hist)):
        g.append(round(G*hist[i]/C*np.exp(-i/T)))
    return g

newhist1=[boltzmaN(i,1) for i in histoT1]
newhist2=[boltzmaN(i,10) for i in histoT10]
print(newhist1)
print(newhist2)
mybin=[0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11]
it=1
for i in newhist1:
        plt.figure()
        plt.title("T=1, Sequence:{}".format(it))
        plt.xlabel("Energy (E)")
        plt.ylabel("Number of states (g)")
        pos = np.arange(len(mybin))
        width = 1.0     # gives histogram aspect to the bar diagram
        ax = plt.axes()
        ax.set_xticks(pos)
        ax.set_xticklabels(mybin)
        plt.bar(pos,i,1)
        plt.savefig("Histogram_T=1,Sequence{}".format(it))
        it+=1
for i in newhist2:
        plt.figure()
        plt.title("T=1, Sequence:{}".format(it))
        plt.xlabel("Energy (E)")
        plt.ylabel("Number of states (g)")
        pos = np.arange(len(mybin))
        width = 1.0     # gives histogram aspect to the bar diagram
        ax = plt.axes()
        ax.set_xticks(pos)
        ax.set_xticklabels(mybin)
        plt.bar(pos,i,1)
        plt.savefig("Histogram_T=10,Sequence{}".format(it))
        it+=1

