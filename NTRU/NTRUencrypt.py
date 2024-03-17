import numpy as np
import sys
from sympy import Poly, symbols
from NTRU.NTRUutil import *

class NTRUencrypt:
    
    def __init__(self, N=503, p=3, q=256, d=18):
        self.N = N 
        self.p = p
        self.q = q

        self.dr = d
        
        self.g = np.zeros((self.N,), dtype=int)
        self.h = np.zeros((self.N,), dtype=int)
        self.r = np.zeros((self.N,), dtype=int)
        self.genr()
        self.m = np.zeros((self.N,), dtype=int)
        self.e = np.zeros((self.N,), dtype=int)
        
        self.I         = np.zeros((self.N+1,), dtype=int)
        self.I[self.N] = -1
        self.I[0]      = 1

        self.readKey = False

        self.Me = None 
        

    def readPub(self,filename="key.pub"):
        with open(filename,"r") as f:
            self.p  = int(f.readline().split(" ")[-1])
            self.q  = int(f.readline().split(" ")[-1])
            self.N  = int(f.readline().split(" ")[-1])
            self.dr = int(f.readline().split(" ")[-1])
            self.h  = np.array(f.readline().split(" ")[3:-1],dtype=int)
        self.I         = np.zeros((self.N+1,), dtype=int)
        self.I[self.N] = -1
        self.I[0]      = 1
        self.genr()
        self.readKey = True


    def genr(self):
        self.r = genRand10(self.N,self.dr,self.dr)
        

    def setM(self,M):
        """
        Set the class message M after performing error checks.
        Before calling this the public key values must have been set (i.e. read)
        NOTE : Message M must be an array describing polynomial coefficients, where the
               polynomial must be degree < N.
        NOTE : The coeffcients must be in [-p/2,p/2].
        NOTE : Message array must be an integer array.
        """
        if self.readKey==False:
            sys.exit("ERROR : Public key not read before setting message")
        if len(M)>self.N:
            sys.exit("ERROR : Message length longer than degree of polynomial ring ideal")
        for i in range(len(M)):
            if M[i]<-self.p/2 or M[i]>self.p/2:
                sys.exit("ERROR : Elements of message must be in [-p/2,p/2]")
        self.m = padArr(M,self.N)

            
    def encrypt(self,m=None):
        """
        Encrypt the message m into the array e
        NOTE : The message m must be set before this routine is called
        """
        if self.readKey == False:
            sys.exit("Error : Not read the public key file, so cannot encrypt")
        if m is not None:
            if len(m)>self.N:
                sys.exit("\n\nERROR: Polynomial message of degree >= N")
            self.m = m
        x = symbols('x')
        self.e = np.array(((((Poly(self.r,x)*Poly(self.h,x)).trunc(self.q)) \
                            + Poly(self.m,x))%Poly(self.I,x)).trunc(self.q).all_coeffs(), dtype=int )
        self.e = padArr(self.e,self.N)
        

    def encryptString(self,M):
        """
        Encrypt the input string M by first converting to binary

        NOTE : The public key must have been read before running this routine
        """

        if self.readKey == False:
            sys.exit("Error : Not read the public key file, so cannot encrypt")
        
        bM = str2bit(M)
        bM = padArr(bM,len(bM)-np.mod(len(bM),self.N)+self.N)
        
        self.Me = ""

        for E in range(len(bM)//self.N):
            self.genr()                               
            self.setM(bM[E*self.N:(E+1)*self.N])      
            self.encrypt()                            
            self.Me = self.Me + arr2str(self.e) + " "

        
