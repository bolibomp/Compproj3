from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import time
#import Boltzman
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
    
def f(x,N):
    return (N-x)/N

start = time.time()
sumE=[]
M=20
states=[1,1,1,14,10]
Eground=[-10,-11,-9,-10,-7]
tempf=[0.5,1,2,10]

S = [[1,1,0,1,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,1,0,1,0,1],
     [0,1,1,1,1,0,1,1,1,0,0,0,0,0,1,1,1,0,1,1,1,0,1,1,1],
     [0,1,1,1,0,0,0,1,0,0,1,1,0,1,1,0,0,0,0,0,1,0,0,1,1],
     [1,0,0,1,0,1,0,1,1,1,0,0,0,0,0,1,0,1,0,1,0,1,1,0,1],
     [1,0,0,0,1,0,1,0,1,0,0,1,0,0,0,1,1,1,0,0,0,1,1,0,0]]
N=30000


 #PHPPHxHxHHHHHPxHxHHPPH=[0,1,0,0,1,x,1,x,1,1,1,1,1,0,x,1,x,1,1,0,0,1]
 #HPPHPHPHHHPPPPPHPHPHPHHPH=[1,0,0,1,0,1,0,1,1,1,0,0,0,0,0,1,0,1,0,1,0,1,1,0,1]
scoreTinf=[0,0,0,0,0]
a=0
for seq in S: #T=infinity
    for n in range(1):
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
        

scoreexpo=[[[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
           [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
           [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
           [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]]


for x in range(0,4):
    a=0
    for T in tempf:
        b=0
        for seq in S:
            for n in range(M):  #T=
                E=0 #[0,1,1,1,2] # x positions for sequence
                Enew=0
                Tnow=T
                
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
                    if Enew==Eground[b]:
                        scoreexpo[x][a][b]+=1
                        break
                    if accept(E,Enew,Tnow)==True:
                        pos=s.new_stat
                        E=Enew
                    Tnow=Tnow*f(x+1,N)
            b+=1
        a+=1
        print(x,T)

finalscoreTinf=np.multiply(scoreTinf,1/M)
finalscoreexpo=np.multiply(scoreexpo,1/M) 

for i in range(len(S)):  
    plt.figure()
    plt.title("Optimization, Sequence: {}".format(i))
    plt.xlabel("Temperature (K)")
    plt.ylabel("Performance number")
    plt.axhline(finalscoreTinf[i], color="red",linestyle="--", label="0 temperature")
    plt.plot(tempf,[x[0][i] for x in finalscoreexpo],"go",label="Annealing x=1")
    plt.plot(tempf,[x[1][i] for x in finalscoreexpo],"bo",label="Annealing x=2")
    plt.plot(tempf,[x[2][i] for x in finalscoreexpo],"ro",label="Annealing x=3")
    plt.plot(tempf,[x[3][i] for x in finalscoreexpo],"yo",label="Annealing x=4")
    plt.legend()
    plt.savefig("Optimization_Sequence{}".format(i))
#for i in range(len(sumE)):
#    if sumE[i]!=0:
#        k+=1
#score1=k/M
#score2=1/(states*M)*np.sum(sumE)
end = time.time() 
print(end-start)
#print(score1)
#print(score2)

        

    # TODO: Implement your metropolis algorithm