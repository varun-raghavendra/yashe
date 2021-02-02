import random
import numpy as np
import math

phid = [1, 0, 1]
q = 29
t = 5
probabilities = []
pi = 3.141592653589793
exp= 2.718281828459045
B = 18
sigma = np.floor(B//6)

def reduce(p):
    global phid
    global q
    p[:] = [x % q for x in p]
    # print(p)
    Q, R = np.polydiv(p, phid)
    return R

def reduce_t(p):
    global phid
    global t
    p[:] = [x % t for x in p]
    # print(p)
    Q, R = np.polydiv(p, phid)
    return R

def recurse(pos, d, myPol, polynomials):
#     global polynomials
    global q
    if pos == d+1:
        # print(myPol)
        newList = [int(x) for x in myPol]
        newList = reduce(newList)
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
        x = exp**(-(np.pi)*(i**2)/sigma**2)
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
    polynomials = []
    genPolynomials(d/2, polynomials)
    return d, polynomials

def inverse(f, polynomials):
    for p in polynomials:
        f_f_inv = (reduce(np.polymul(p, f))).tolist()
        if f_f_inv == [1]:
#             print(p)
            return p
    return [-1]

def KeyGen(d, polynomials):
    global phid
    global q
    global t
    f_prime = ChiKey().tolist()
    g = ChiKey().tolist()
    f = [x*t for x in f_prime]
    f[-1] += 1
    reduce(f)
    
    f_inv = []
    while True:
        X = (inverse(f, polynomials))
        f_inv = [int(x) for x in X]
        if f_inv != [-1]:
            break
    
    h = reduce(np.polymul(g, f_inv))
    print(h)
    h_prime = [x*t for x in h]
    h = reduce(h_prime)
    # print(h)
    
    return f_prime, g, f_inv, reduce(h.tolist()), reduce(f)

def Encrypt(h, msg):
    global q
    global t
    e = ChiErr(2)
    s = ChiErr(2)
    # print("e = ", e)
    # print("s = ", s)
    
    delta = math.floor(q/t)
    c = [delta*x for x in msg]
    c = np.polyadd(c, e)
    H = reduce(np.polymul(h, s))
    c = np.polyadd(c, H)
    return reduce(c)

def Decrypt(f, c):
    # print(np.polymul(f, c))
    f_c = np.polymul(f, c)
    Q1, R1 = np.polydiv(f_c, phid)
    print(R1)
    M = reduce(np.polymul(f, c)).tolist()
    # print(M)
    M[:] = [round(x*t/q) for x in M]

    return reduce_t(M)
    
d, polynomials = ParamsGen()
genProbabilities()
# print(d, phid, q, t)

f_prime, g, f_inv, h, f = KeyGen(d, polynomials)
f = reduce(f)
print(f_prime, g, f_inv, h, f)

print("Enter message (a, b) to encrypt")
s = input()
msg = s.split()
msg[:] = [int(x) for x in msg]
print("The message entered is:")
print(msg)

c = Encrypt(h, reduce_t(msg))
print("Cipher text obtained is")
print(c)

final_msg = Decrypt(f, c)
print("Final message obtained is")
print(final_msg)

# e = ChiErr