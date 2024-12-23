" Created by  Guillermo Romero University of Santiago of Chile USACH, guillermo.romero@usach.cl"
 

from scipy import *
from scipy import sparse
from pylab import *
import scipy
import pickle
import matplotlib.pyplot as plt
import numpy 
#import time
import cmath
import math
from matplotlib.font_manager import FontProperties
import matplotlib.patches as patches
from csv import *
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import ticker, cm
from scipy.linalg import expm, sinm, cosm, logm
from numpy import linalg as LA
from tqdm import tqdm


def BH_exact_dynamics(sites, h, J, N):

    "Defining local operators"
    "Initialize local operators"
    L = sites
    #if L >= 6:
     #   warnings.warn("Large L: Exact diagonalization might take a long time!")
    Energy = []
    id = eye(2)
    ID = eye(2**sites)
    sx = np.array([[0., 1.],[1., 0.]])
    sy = np.array([[0., -1j],[1j, 0.]])
    sz = np.array([[1., 0.],[0., -1.]])
    #print(am,ap)
    sx_list = []  # sx_list[i] = kron([id, id, ..., id, sx, id, .... id])
    sz_list = []
    dim = 2**L
    for i_site in range(L):
        sx_ops = [id] * L
        sz_ops = [id] * L
        sx_ops[i_site] = sx
        sz_ops[i_site] = sz
        Sx = sx_ops[0]
        Sz = sz_ops[0]
        for j in range(1, L):
            Sx = numpy.kron(Sx, sx_ops[j])
            Sz = numpy.kron(Sz, sz_ops[j])
        sx_list.append(Sx)
        sz_list.append(Sz)
              
    H_xx = numpy.zeros((dim, dim))
    H_z = numpy.zeros((dim, dim))
    H_x = numpy.zeros((dim, dim))

    for i in range(L - 1):
        H_xx = H_xx + sx_list[i].dot(sx_list[(i + 1) % L]) 
        #H_xx = H_xx + ap_list[i] * am_list[(i + 1) % L] + am_list[i] * ap_list[(i + 1) % L]
    for i in range(L):
        H_z = H_z + h * sz_list[i]
        H_x = H_x + sx_list[i]
    H = J * H_xx  + H_z
    
    D, V = LA.eigh(H)
    D = sorted(D, reverse=False) # Sorted descending
    idx = numpy.argsort(D)
    V = V[:,idx]
  

    "Oscillating driving and dynamical dimer formation"
    "Computing the ground state with filling factor n = 1"
    "Quantum states"
   
    Fock = scipy.sparse.eye(2).toarray()
    
    p1 = Fock[:,0]
    p0 = Fock[:,1]

    
    P0 = numpy.kron(numpy.kron(p0,p0),p0)
    P1 = numpy.kron(numpy.kron(p1,p1),p0)
    P2 = numpy.kron(numpy.kron(p1,p0),p1)
    P3 = numpy.kron(numpy.kron(p0,p1),p1)
    
      
    "Stroboscopic dynamics"
    "Integer"
    Omega0 = 2*h
    omega1 = Omega0 
    omega2 = 2*Omega0
   
    T = 2*pi/(Omega0)
    t_array1 = np.array([(omega2 if i % 2 == 0 else omega1) for i in range(L-1)])

  
    scale = 1/h
    steps = 2000
    dt1 = T/steps
    psi0 = p0
    print(dt1 * scale)
    for i in range(L-1):
        psi0 = kron(psi0,p0)
        
    #psi0 = psi0
    tspan1 = numpy.linspace(0.0,T,steps)
    U0 = ID
    
    
    for z in range(len(tspan1)):
        H1 = numpy.zeros((dim, dim))
        for ii in range(L-1):
            H1= H1 + J * numpy.cos(t_array1[ii]*tspan1[z]) * sx_list[ii].dot(sx_list[ii+1])
        HBH = H_z + H1
        Udt = expm(-1j * HBH * dt1 * scale)
        U1 = numpy.matmul(Udt,U0)
        U0 = U1
    
    psi0 =p1
   
    for i in range(L-1):
        psi0 = kron(psi0,p1)
        
    psi0 = psi0
    time = numpy.arange(1.,N+1)  

    VnE = numpy.zeros((len(time)))
    SZT = numpy.zeros((len(time),L))
    SXT = numpy.zeros((len(time),L))
    SZSZ = numpy.zeros((len(time),L-1))

    
    for p in tqdm(range(len(time))):
        psiT = U1.dot(psi0)
                   
        
        for ii in range(L):
            SZT[p,ii] = psiT.conjugate().transpose().dot(sz_list[ii].dot(psiT))
            SXT[p,ii] = psiT.conjugate().transpose().dot(sx_list[ii].dot(psiT))
        for ii in range(L-1):
            SZSZ[p,ii] = psiT.conjugate().transpose().dot(sz_list[ii].dot(sz_list[ii+1]).dot(psiT))-psiT.conjugate().transpose().dot(sz_list[ii].dot(psiT))*psiT.conjugate().transpose().dot(sz_list[ii+1].dot(psiT))
               
        psi0 = psiT


   
    "COLORS AND MARKERS"
    
    colors = ('tab:blue','tab:orange','tab:red','tab:green','tab:green','k','tab:olive','tab:brown','tab:orange','#99CCFF','#004C99','#FFCC99')

    "Blues"
    #colors = ('#CCE5FF','#99CCFF','#66B2FF','#3399FF','#0080FF','#0066CC','#004C99','#FF9999','#CC0000')
    markers = ('+','o','v','^','<','>','s','h','x','None','*')


    font = {'family' : 'serif',
            'serif': 'Times New Roman', 
            'weight' : 'regular',
            'style' : 'italic',
            
            'size'   : 20}

    plt.rcParams['axes.linewidth'] = 1.5 #set the value globally
    matplotlib.rc('font', **font)


    fig = plt.figure(figsize=(10,6))
    ax1 = fig.add_subplot(1,1,1)
    plt.plot(time,SZT[:,0],'o',time,SZT[:,1],'s',time,SZT[:,2],'x')
   
    fig = plt.figure(figsize=(10,6))
    ax2= fig.add_subplot(1,1,1)
   
    plt.plot(time,SZSZ[:,0],'s',time,SZSZ[:,1],'o')
    plt.show()
  
    

   
if __name__ == "__main__":
    sites = 3
    h = 1. # Transverse field
    J = h/10 # Spin-spin exchange
    N = 50 # Number of periods for the stroboscopiv dynamics
    BH_exact_dynamics(sites=sites, h = h, J = J, N=N)
 
  
