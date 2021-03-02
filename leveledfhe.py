from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_gcdex
import random
import numpy as np
import math

phid = [1, 0, 0, 0, 0, 0, 0, 0, 1]
q = 1040101
t = 50
w = 4
probabilities = []
pi = 3.141592653589793
exp= 2.718281828459045
B = 18
sigma = np.floor(B/6)

def reduce(p, mod, centered=True, divide=True):
    # Put the given polynomial in the range [-mod/2, mod/2) if centered,
    # and to [0,mod) if not centered.
    global phid

    R = [x % mod for x in p]
    # print(p)
    if divide:
        Q, R = np.polydiv(p, phid)
        R = [x % mod for x in R]

    if centered:
        for i in range(len(R)):
            if R[i]>mod//2:
                R[i]-=mod

    return R

def genProbabilities():
    # Make probability space assuming B is our upper limit for Chi
    global probabilities
    global B
    sum = 0
    for i in range (-1*B+1, B):
        x = np.e**(-(np.pi)*(i**2)/sigma**2)
        probabilities.append(x)
        sum += x
    probabilities = [x / sum for x in probabilities]

def sample(B):
    # Returns a random sample from -B+1 to B, according to the
    # probability space
    return np.random.choice(np.arange(-1*B+1, B), p=probabilities)

def ChiErr(n, B_err = B):
    # Makes a small coefficient error polynomial from the Chi function
    poly = []
    for i in range (0, n):
        ele = sample(B_err)
        poly.append(ele)
    return poly

def ChiKey(n, B_key = B):
    # Makes a small coefficient key polynomial from the Chi function
    poly = []
    for i in range (0, n):
        ele = sample(B_key)
        poly.append(ele)
    return poly

def inverse(f):
    # Fast inverse method
    global phid
    global q
    p = ZZ.map(f)
    mod = ZZ.map(phid)
    s, t, g = gf_gcdex(p, mod, q, ZZ)
    if len(g) == 1 and g[0] == 1:
        return s
    else:
        return [-1]

def ParamsGen(Lambda=987654321):
    """
    Given security parameter lambda, initialise d, q, t, ChiKey,
    ChiErr, and w, where w is an integer > 1
    """
    d = 16
    w = sample(B)

    return d, w

def KeyGen():
    d, w = ParamsGen()

q = 7
w = 2

def BitDecomp(x):
    Lwq = int(np.floor(np.log(q)/np.log(w)))+2
    # Lwq = 3
    res = []
    x_bar = x
    for i in range(1, Lwq):
        x_bar = [xi* w for xi in x_bar]
        x_bar = reduce(x_bar, q, centered=True, divide=False)
        print(type(x_bar))
        res.append(x_bar)

    print(res)


def pow_of_2(x):
    #global lwq
    global w
    lwq =  math.floor(math.log(q) / math.log(w)) + 2

    z = x
    y = []
    for i in range(lwq):
        z = [j*(w**i) for j in x]
        #print(i,z,reduce(z))
        z=reduce(z)
        if z==[0.0]:
            z=z*len(x)
        y.append(z)
    return y


BitDecomp([1, 2])
