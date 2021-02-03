from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_gcdex
import random
import numpy as np
import math

phid = [1, 0, 1]
q = 2003
t = 7
probabilities = []
pi = 3.141592653589793
exp= 2.718281828459045
B = 18
sigma = np.floor(B/6)

def reduce(p):
    global phid
    global q
    p[:] = [x % q for x in p]
    # print(p)
    Q, R = np.polydiv(p, phid)
    R[:] = [x % q for x in R]
    for i in range(len(R)):
      if R[i]>q//2:
        R[i]-=q
    return R

def reduce_t(p):
    global phid
    global t
    p[:] = [x % t for x in p]
    # print(p)
    Q, R = np.polydiv(p, phid)
    R[:] = [x % t for x in R]
    for i in range(len(R)):
      if R[i]>t//2:
        R[i]-=t
    return R

def reduceMod(p):
    global phid
    global q
    p[:] = [x % q for x in p]
    # print(p)
    Q, R = np.polydiv(p, phid)
    R[:] = [x % q for x in R]
    
    return R

def reduceMod_t(p):
    global phid
    global t
    p[:] = [x % t for x in p]
    # print(p)
    Q, R = np.polydiv(p, phid)
    R[:] = [x % t for x in R]
    
    return R

def recurse(pos, d, myPol, polynomials):
#     global polynomials
    global q
    if pos == d+1:
        # print(myPol)
        newList = [int(x) for x in myPol]
        newList = reduceMod(newList)
        polynomials.append(newList)
        return
    for i in range (0, q):
        myPol.append(i)
        recurse(pos+1, d, myPol, polynomials)
        myPol.pop()

def genPolynomials(d, polynomials):
    global q
    myPol = []
    recurse(0, d, myPol, polynomials)

def genProbabilities():
    global probabilities
    global B
    sum = 0
    for i in range (-1*B+1, B):
        x = np.e**(-(np.pi)*(i**2)/sigma**2)
        probabilities.append(x)
        sum += x
    probabilities[:] = [x / sum for x in probabilities]

def sample():
    return np.random.choice(np.arange(-1*B+1, B), p=probabilities)

def ChiErr(n):
    poly = []
    for i in range (0, n):
        ele = sample()
        poly.append(ele)
    return poly

def ChiKey():
    global q
    coeff = np.array([random.randint(1,q-1) for i in range((d>>1)+1)])
    return coeff

def ParamsGen():
    d = 4
    return d

def inverse(f):
    global phid
    global q
    p = ZZ.map(f)
    mod = ZZ.map(phid)
    s, t, g = gf_gcdex(p, mod, q, ZZ)
    if len(g) == 1 and g[0] == 1: 
        return s
    else:
        return [-1]

# def inverse(f, polynomials):
#     for p in polynomials:
#         f_f_inv = (reduceMod(np.polymul(p, f))).tolist()
#         if f_f_inv == [1]:
# #             print(p)
#             return p
#     return [-1]

def KeyGen(d):
    global phid
    global q
    global t
    f_prime = ChiErr(d)
    g = ChiErr(d)
    f = [x*t for x in f_prime]
    f[-1] += 1
    reduce(f)

    print("f_prime = ", f_prime)
    print("g = ", g)
    print("f = ", f)
    f_inv = []
    while True:
        X = (inverse(f))
        f_inv = [int(x) for x in X]
        if f_inv != [-1]:
            break

    h = reduce(np.polymul(g, f_inv))
    h_prime = [x*t for x in h]
    h = reduce(h_prime)

    return f_prime, g, f_inv, reduce(h.tolist()), reduce(f)

def Encrypt(h, msg):
    global q
    global t
    e = reduce(ChiErr(2))
    s = reduce(ChiErr(2))
    # print("e = ", e)
    # print("s = ", s)

    delta = math.floor(q/t)
    c = [delta*x for x in msg]
    c = np.polyadd(c, e)
    H = reduce(np.polymul(h, s))
    c = np.polyadd(c, H)
    return reduce(c)

def Decrypt(f, c):
    M = reduce(np.polymul(f, c)).tolist()
    # print(M)
    M[:] = [round(x*t/q) for x in M]

    return reduce_t(M)

def HomomorphicAddition(c1, c2):
    A = np.polyadd(c1, c2)
    return reduce(A).tolist()

d = ParamsGen()
genProbabilities()
# print(d, phid, q, t)

f_prime, g, f_inv, h, f = KeyGen(d)
# f_prime = [0, 1, -1, 1]
# g = [-1, -1, 0, 0]
# f = [-4, 1]
# h = [-40, -36]
# f_inv = [12, 53]
print("f_prime =",f_prime,"g =", g,"f_inv =", f_inv,"h =", h,"f =", f)
f = reduce(f)
print("f_prime =",f_prime,"g =", g,"f_inv =", f_inv,"h =", h,"f =", f)

print("Enter message (a, b) to encrypt")
s = input()
msg = s.split()
msg[:] = [int(x) for x in msg]
print("The message entered is:")
print(msg)

c = Encrypt(h, reduce_t(msg))
print("Cipher text obtained is")
print(c)

delta = math.floor(q/t)
del_msg = [delta*x for x in msg]
print("inherent noise in encryption is:")
print(reduce(np.polysub(np.polymul(f, c), del_msg)))

final_msg = Decrypt(f, c)
print("Final message obtained is")
print(final_msg)

# Homomorphic Addition
print("Enter two messages (a1, b1) and (a2, b2) to encrypt")
s1 = input()
s2 = input()
msg1 = s1.split()
msg2 = s2.split()
msg1[:] = [int(x) for x in msg1]
msg2[:] = [int(x) for x in msg2]
print("The messages entered are:")
print("message1: ", msg1)
print("message2: ", msg2)

c1 = Encrypt(h, reduce_t(msg1))
c2 = Encrypt(h, reduce_t(msg2))

final_msg = Decrypt(f, HomomorphicAddition(c1, c2))
print(final_msg)
# e = ChiErr
