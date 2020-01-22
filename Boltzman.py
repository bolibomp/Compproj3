from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import time

from scipy.spatial.distance import pdist, squareform



class Switcher(object):
    """
    a class to switch cases and calculate the hh contacts in each SAW sequence
    """
    
    def __init__(self,seq, xpos, ypos):

        self.seq = seq      # HP list which contains of  0=P & 1=H 
        self.xpos = xpos    # x_position of elements
        self.ypos = ypos    # y_position of elements
        self.n = len(seq)

        # any_in will be used to ckeck SAW (self-avoiding walk)
        self.any_in = lambda a, b: any(i in b for i in a)
        # calculating distance
        self.dist = lambda d : np.sqrt((d ** 2).sum(axis=1))

        
    def pivot(self):
        
        k = np.random.randint(1, self.n - 1) # random pivot point 

        # set the origin to (0,0) for chunk[k+1:]
        x_tmp = [a - self.xpos[k] for a in  self.xpos[k+1:]]
        y_tmp = [b - self.ypos[k] for b in  self.ypos[k+1:]]
        n_tmp = len(x_tmp)
         
        
        g = np.random.randint(0, 7) # random rotation cases
        method_name = 'case_' + str(g)
        method = getattr(self, method_name, lambda: "nothing")
        
        # transform the chunk[k+1:]
        x_rot,y_rot = method(x_tmp, y_tmp, n_tmp )
        


        # move back to the origin
        x_new = self.xpos[:k+1] + [a + self.xpos[k] for a in  x_rot]
        y_new = self.ypos[:k+1] + [b + self.ypos[k] for b in  y_rot]
        
        new_stat = list(zip(x_new, y_new))

             
        
        # check if the new state is a SAW      
        if self.any_in(new_stat[:k+1], new_stat[k+1:]):           
            # recursive pivot - until it is saw 
            # incase of keeping the old state comment the next line
            return self.pivot()   
        else:
            # if it is SAW calculate HH contacts
            self.new_stat = new_stat
            return('HH contacts  ',self.hh_contact())
                    

    
    def hh_contact(self):
        #calculate HH contacts
        seq = np.array(self.seq)
        n_stat = np.array(self.new_stat)
        dist =  squareform( pdist(n_stat, metric='euclidean'))
        arr =np.triu (dist, 2)
        index = np.argwhere(arr==1)
        E = 0
        for c in index:
            if seq[c].all()==1:
                E-=1        
        return E

      
    
    @staticmethod
    def case_0( x, y, n):        #+90 rotation
        for i in range(n):
            tmp = x[i]
            x[i] = - y[i]
            y[i] = tmp
        return  x , y
    
    @staticmethod
    def case_1(x, y, n):        #-90 rotation
        for i in range(n):
            tmp = x[i]
            x[i] =  y[i]
            y[i] = -tmp 
        return x ,  y
    
    @staticmethod   
    def case_2(x, y, n):        #180 rotation
        for i in range(n):
            x[i] = - x[i]
            y[i] = - y[i]
        return  x , y
        
    @staticmethod    
    def case_3(x, y, n):        # mirror(1,0)
        for i in range(n):
            y[i] = - y[i]
        return  x , y

    @staticmethod    
    def case_4(x, y, n):        # mirror(0,1)
        for i in range(n):
            x[i] = - x[i]
        return   x ,  y

    @staticmethod
    def case_5(x, y, n):        # mirror(1,1)
        for i in range(n):
            tmp = x[i]
            x[i] =  y[i]
            y[i] =  tmp
        return   x ,  y
    
    @staticmethod        
    def case_6(x, y, n):        # mirror(1,-1)
        for i in range(n):
            tmp = x[i]
            x[i] = - y[i]
            y[i] = - tmp
        return   x , y


def boltzmaN(hist,T):
    G=5768299665 #total number of states
    u=[]
    g=[]
    for i in range(len(hist)):
        u.append(hist[i]*np.exp(i/T))
    C=np.sum(u)    
    for i in range(len(hist)):
        g.append(G*hist[i]/C*np.exp(i/T))
    return g

def accept(E1,E2,T): 
    if E2<=E1:
        return True
    kb=1
    p=np.exp(-(E2-E1)/(kb*T))
    r=np.random.rand()
    if r<p:
        return True
    else:
        return False

#---------------------------------------------- Task 1 / Producing histograms
M=1 #Runnung the program M times
N=3000  #5000000 #Folding N times
states=[1,1,1,14,10]
Eground=[-10,-11,-9,-10,-7]
temp=[1,10]
S = [[1,1,0,1,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,1,0,1,0,1], #sequences, P=0 H=1
     [0,1,1,1,1,0,1,1,1,0,0,0,0,0,1,1,1,0,1,1,1,0,1,1,1],
     [0,1,1,1,0,0,0,1,0,0,1,1,0,1,1,0,0,0,0,0,1,0,0,1,1],
     [1,0,0,1,0,1,0,1,1,1,0,0,0,0,0,1,0,1,0,1,0,1,1,0,1],
     [1,0,0,0,1,0,1,0,1,0,0,1,0,0,0,1,1,1,0,0,0,1,1,0,0]]

for T in temp: 
    Etot=[]
    l=0
    for seq in S:
        Esave=[0] 
        E=0
        pos=[]
        for i in range(len(seq)):
            pos.append((0,i))
        xpos=[0]*len(seq)
        ypos=[0]*len(seq)
        for k in range(N):
            for i in range(len(seq)):
                xpos[i]=pos[i][0]
                ypos[i]=pos[i][1]
            s=Switcher(seq,xpos,ypos)
            s.pivot()  
            Enew=s.hh_contact()
            if accept(E,Enew,T)==True:
                pos=s.new_stat   
                E=Enew
            Esave.append(Enew)
        Etot.append(Esave)
        print(l)
        l+=1
    
    histtot=[]
    for i in range(len(Etot)):
        hist=[0]*12
        for j in range(N):
            E=np.abs(Etot[i][j])
            hist[E]+=1
        histtot.append(hist)
    
    
    
    newhist=[boltzmaN(i,T) for i in histtot]
    mybin=[0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11]
    it=1
    for i in newhist: #plotting
        plt.figure()
        plt.title("T={}, Sequence:{}".format(T,it))
        plt.xlabel("Energy (E)")
        plt.ylabel("Number of states (g)")
        pos = np.arange(len(mybin))
        width = 1.0
        ax = plt.axes()
        ax.set_xticks(pos)
        ax.set_xticklabels(mybin)
        plt.bar(pos,i,1)
        plt.savefig("Histogram_T={},Sequence:{}.png".format(T,it))
        it+=1
#------------------------------------- Task 3 / Optimization
sumE=[]
M=20 #Running the program M times
N=30000 #Folding N times
states=[1,1,1,14,10]
Eground=[-10,-11,-9,-10,-7]
temp=[0.1,1,10]
S = [[1,1,0,1,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,1,0,1,0,1],
     [0,1,1,1,1,0,1,1,1,0,0,0,0,0,1,1,1,0,1,1,1,0,1,1,1],
     [0,1,1,1,0,0,0,1,0,0,1,1,0,1,1,0,0,0,0,0,1,0,0,1,1],
     [1,0,0,1,0,1,0,1,1,1,0,0,0,0,0,1,0,1,0,1,0,1,1,0,1],
     [1,0,0,0,1,0,1,0,1,0,0,1,0,0,0,1,1,1,0,0,0,1,1,0,0]]


scoreTinf=[0,0,0,0,0]
a=0
for seq in S: #T=0
    for n in range(M):
        E=0
        pos=[]
        for i in range(len(seq)):
            pos.append((0,i))
        xpos=[0]*len(seq)
        ypos=[0]*len(seq)
        
        for k in range(N):
            for i in range(len(seq)):
                xpos[i]=pos[i][0]
                ypos[i]=pos[i][1]
            s=Switcher(seq,xpos,ypos)
            s.pivot()
            Enew=s.hh_contact()
            if Enew==Eground[a]:
                scoreTinf[a]+=1
                break
            if Enew>=E:
                pos=s.new_stat
                E=Enew
    a+=1

scoreexpo=[[0]*len(S) for _ in range(len(temp))]
a=0
for T in temp: #Exponential annealing
    b=0
    for seq in S:
        for n in range(M):
            E=0 
            Enew=0
            pos=[]
            for i in range(len(seq)):
                pos.append((0,i))
            xpos=[0]*len(seq)
            ypos=[0]*len(seq)
            Tnow=T
            
            for k in range(N):
                Tred=Tnow/N
                for i in range(len(seq)):
                    xpos[i]=pos[i][0]
                    ypos[i]=pos[i][1]
                s=Switcher(seq,xpos,ypos)
                s.pivot()
                Enew=s.hh_contact()
                if Enew==Eground[b]:
                    scoreexpo[a][b]+=1
                    break
                if accept(E,Enew,Tnow)==True:
                    pos=s.new_stat
                    E=Enew
                Tnow-=Tred
        b+=1
    a+=1

finalscoreTinf=np.multiply(scoreTinf,1/M)
finalscoreexpo=np.multiply(scoreexpo,1/M) 

for i in range(len(S)):  #Plotting
    plt.figure()
    plt.title("Optimization, Sequence: {}".format(i))
    plt.xlabel("Temperature (K)")
    plt.ylabel("Performance number")
    plt.axhline(finalscoreTinf[i], color="red",linestyle="--", label="0 temperature")
    plt.semilogy(temp,[x[i] for x in finalscoreexpo],"go",label="Exponential Annealing")
    plt.legend()
    plt.savefig("Optimization_Sequence{}".format(i))
    
#---------------------------------------- Task 4 / Finding designing sequences

poslist=[[] for _ in range(16)]
deglist=[[] for _ in range(16)]
M=10 #Running the program M times
N=30000 #Folding N times
k=0
for a in range(2): #producing all binary combinations of a,b,c,d
    for b in range(2):
        for c in range(2):
            for d in range(2):
                for n in range(M):
                    deg=[]
                    seq = [0,1,0,0,1,a,1,b,1,1,1,1,1,0,c,1,d,1,1,0,0,1]
                    E=0
                    pos=[]
                    for i in range(len(seq)):
                        pos.append((0,i))
                    xpos=[0]*len(seq)
                    ypos=[0]*len(seq)
                    T=1
                    Tred=T/N 
                    for n in range(N):
                        for i in range(len(seq)):
                            xpos[i]=pos[i][0]
                            ypos[i]=pos[i][1]
                        s=Switcher(seq,xpos,ypos)
                        s.pivot()
                        Enew=s.hh_contact()
                        if Enew==-11:
                            poslist[k].append(s.new_stat)
                        if accept(E,Enew,T)==True:
                            pos=s.new_stat
                            E=Enew
                        T=T-Tred
                k+=1
for j in range(16):
    for i in poslist[j]:
        if i not in deglist[j]:
            deglist[j].append(i)

degstates=[len(x) for x in deglist] # result
