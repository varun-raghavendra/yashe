from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_gcdex
import random
import numpy as np
import math

"""
Things we can polish:

1. Take variables out of global scope,
2. Shrink reduce, reduce_t to one function,
3. Fix a ratio between B and q
4. Make function to derive PhiD for a given d
5. Make a demo function that formats output neatly, for final demo.

"""

phid = [1, 0, 0, 0, 0, 0, 0, 0, 1]
q = 1040101
t = 50
probabilities = []
pi = 3.141592653589793
exp= 2.718281828459045
B = 18
sigma = np.floor(B/6)

# def reduce(p):
#     # Put the given polynomial in the range [-q/2, q/2)
#     global phid
#     global q
#     p = [x % q for x in p]
#     # print(p)
#     Q, R = np.polydiv(p, phid)
#     R = [x % q for x in R]
#     for i in range(len(R)):
#       if R[i]>q//2:
#         R[i]-=q
#     return R
#
# def reduce_t(p):
#     # Put the given polynomial in the range [-t/2, t/2)
#     global phid
#     global t
#     p = [x % t for x in p]
#     # print(p)
#     Q, R = np.polydiv(p, phid)
#     R = [x % t for x in R]
#     for i in range(len(R)):
#       if R[i]>t//2:
#         R[i]-=t
#     return R
#
# def reduceMod(p):
#     # Put the given polynomial in the range [-0, q)
#     global phid
#     global q
#     p = [x % q for x in p]
#     # print(p)
#     Q, R = np.polydiv(p, phid)
#     R = [x % q for x in R]
#
#     return R
#
# def reduceMod_t(p):
#     # Put the given polynomial in the range [-0, t)
#     global phid
#     global t
#     p = [x % t for x in p]
#     # print(p)
#     Q, R = np.polydiv(p, phid)
#     R = [x % t for x in R]
#
#     return R

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


# def recurse(pos, d, myPol, polynomials):
#     # Deprecated function to make every polynomial in Z[X]
#
#     global q
#     if pos == d+1:
#         # print(myPol)
#         newList = [int(x) for x in myPol]
#         newList = reduceMod(newList)
#         polynomials.append(newList)
#         return
#     for i in range (0, q):
#         myPol.append(i)
#         recurse(pos+1, d, myPol, polynomials)
#         myPol.pop()

# def genPolynomials(d, polynomials):
#     # Deprecated Initialiser function to call recurse and
#     # set up polynomial space.
#     global q
#     myPol = []
#     recurse(0, d, myPol, polynomials)

# def inverse(f, polynomials):
#     # Deprecated slow inverse
#     for p in polynomials:
#         f_f_inv = (reduceMod(np.polymul(p, f))).tolist()
#         if f_f_inv == [1]:
# #             print(p)
#             return p
#     return [-1]

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

def ParamsGen():
    # initialiser function.
    d = 16
    return d

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

def KeyGen(d):
    global phid
    global q
    global t
    f_prime = ChiErr(d)
    g = ChiErr(d)
    f = [x*t for x in f_prime]
    f[-1] += 1
    reduce(f,q)

    # print("f_prime = ", f_prime)
    # print("g = ", g)
    # print("f = ", f)
    f_inv = []
    while True:
        try:
            X = (inverse(f))
        except:
            f_prime = ChiErr(d)
            f = [x*t for x in f_prime]
            f[-1] += 1
            reduce(f,q)
            continue
        f_inv = [int(x) for x in X]
        if f_inv != [-1]:
            break

    h = reduce(np.polymul(g, f_inv),q)
    h_prime = [x*t for x in h]
    h = reduce(h_prime,q)

    return f_prime, g, f_inv, reduce(h,q), reduce(f,q)

def Encrypt(h, msg):
    global q
    global t
    e = reduce(ChiErr(9),q) #This shouldn't be 2, it should be some
    s = reduce(ChiKey(9),q) #dependent of d, no?


    print("e = ", e)
    # print("s = ", s)

    delta = math.floor(q/t)
    c = [delta*x for x in msg]
    c = np.polyadd(c, e)
    H = reduce(np.polymul(h, s),q)
    c = np.polyadd(c, H)
    return reduce(c,q)

def Decrypt(f, c):
    global t
    M = reduce(np.polymul(f, c),q)
    # print(M)
    M = [round(x*t/q) for x in M]

    return reduce(M,t)

def HomomorphicAddition(c1, c2):
    global q
    A = np.polyadd(c1, c2)
    return reduce(A,q)

d = ParamsGen()
genProbabilities()
# print(d, phid, q, t)

f_prime, g, f_inv, h, f = KeyGen(d)
# f_prime = [0, 1, -1, 1]
# g = [-1, -1, 0, 0]
# f = [-4, 1]
# h = [-40, -36]
# f_inv = [12, 53]
# print("f_prime =",f_prime,"g =", g,"f_inv =", f_inv,"h =", h,"f =", f)
f = reduce(f,q)
# print("f_prime =",f_prime,"g =", g,"f_inv =", f_inv,"h =", h,"f =", f)

# print("Enter message (a, b) to encrypt")
# s = input()
# msg = s.split()
# msg = [int(x) for x in msg]
# # print("The message entered is:")
# print(msg)

# c = Encrypt(h, reduce_t(msg))
# print("Cipherd text obtained is")
# print(c)

# delta = math.floor(q/t)
# del_msg = [delta*x for x in msg]
# print("inherent noise in encryption is:")
# print(reduce(np.polysub(np.polymul(f, c), del_msg)))

# final_msg = Decrypt(f, c)
# print("Final message obtained is")
# print(reduceMod_t(final_msg))

# Homomorphic Addition
print("Enter two messages (a1, b1) and (a2, b2) to encrypt")
s1 = input()
s2 = input()
msg1 = s1.split()
msg2 = s2.split()
msg1 = [int(x) for x in msg1]
msg2 = [int(x) for x in msg2]
print("The messages entered are:")
print("message1: ", msg1)
print("message2: ", msg2)

c1 = Encrypt(h, reduce(msg1,t))
c2 = Encrypt(h, reduce(msg2,t))

final_msg = Decrypt(f, HomomorphicAddition(c1, c2))
print("First cipher text is: ", c1)
print("Second cipher text is: ", c2)
print("Homomorphic Addition is ", reduce(final_msg,t,centered=False))
# e = ChiErr
