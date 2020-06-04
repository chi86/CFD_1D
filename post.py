#! python3

import copy,os,math,sys
sys.path.append(os.environ['PYTHONPATH']+"/")

import numpy as np
import matplotlib.pyplot as plt

from cell import Cell
from mesh import Mesh
from fluidProps import FluidPopsRef
import rans
import dns
from postProcess import PostProcess

import warnings
warnings.simplefilter("ignore")

def main():
    ## DNS const data
    filenameDNS=os.environ['DNS']+"/re360pr20/"
    DNS_C=dns.readDNSdataConst(filenameDNS)
    
    filenameDNS=os.environ['DNS']+"/Var_Re360Pr20/"
    DNS_V=dns.readDNSdataVar(filenameDNS)

    CFD1D=readCFD1D()
    

    PostProcess([CFD1D,DNS_C,DNS_V],["w","chi","k","epsilon","Zeta","psi","muT","aT","tauL","tauT","qL","qT"])




def readCFD1D():
    data = np.loadtxt("data.dat", skiprows=1)
    imax=len(data)
    CFD1D=Mesh(imax)

    CFD1D.typ="data"

    for idx,val in enumerate(data[:]):
        CFD1D.cells[idx+1].ru=val[0]
        CFD1D.cells[idx+1].rp=val[1]
        CFD1D.cells[idx+1].dr=val[2]
        CFD1D.cells[idx+1].yp=val[3]

        CFD1D.cells[idx+1].w=val[4]
        CFD1D.cells[idx+1].theta=val[5]
        CFD1D.cells[idx+1].chi=val[6]
        
        CFD1D.cells[idx+1].k=val[7]
        CFD1D.cells[idx+1].epsilon=val[8]
        CFD1D.cells[idx+1].Zeta=val[9]
        CFD1D.cells[idx+1].psi=val[10]
        CFD1D.cells[idx+1].V2=val[19]
        CFD1D.cells[idx+1].muT=val[11]
        
        CFD1D.cells[idx+1].L=val[12]
        CFD1D.cells[idx+1].T=val[13]
        
        CFD1D.cells[idx+1].aT=val[14]
        
        CFD1D.cells[idx+1].tauL=val[15]
        CFD1D.cells[idx+1].tauT=val[16]
        
        CFD1D.cells[idx+1].qL=val[17]
        CFD1D.cells[idx+1].qT=val[18]

    return CFD1D


if __name__ == '__main__':

    main()
