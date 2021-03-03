from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_gcdex
import random
import numpy as np
import math
import itertools

phid = [1, 0, 0, 0, 0, 0, 0, 0, 1]
d = 16
q = 1040101
t = 50
w = 4
n = 9
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


def Basic_KeyGen(d):
    global phid
    global q
    global t
    global B

    f_prime = ChiErr(d)
    g = ChiErr(d)
    f = [x*t for x in f_prime]
    f[-1] += 1
    reduce(f, q, centered=True)

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
            reduce(f, q, centered=True)
            continue
        f_inv = [int(x) for x in X]
        if f_inv != [-1]:
            break

    h = reduce(np.polymul(g, f_inv), q, centered=True)
    h_prime = [x*t for x in h]
    h = reduce(h_prime, q, centered=True)

    return f_prime, g, f_inv, reduce(h, q, centered=True), reduce(f, q, centered=True)


def split_to_bits(x, w, Lwq):
    ans = []
    temp = []

    base = w
    temp = np.zeros(Lwq)

    if(x<0):
        # temp[Lwq-1]+=1
        x = pow(w, Lwq)+x

    idx = 0

    while(x):
        t = x%w
        temp[idx]=t
        idx+=1
        x = x//w

    # print(temp)
    return temp

# split_to_bits(-1, 2, 6)

def BitDecomp(x):
    Lwq = int(np.floor(np.log(q)/np.log(w)))+2
    # Lwq = 3
    # res = np.array(0)
    res = []

    for x_coeff in x:
        temp = split_to_bits(x_coeff, w, Lwq)
        res.append(temp)


    # print(res)
    res = np.asarray(res)
    # The res here, is the coefficients themselves as an Lwq array.
    # We won't be adding directly to these, but we'll be adding to the
    # columns of the array. So it makes sense to split it into that
    # form and send.

    # print(res)
    ans = []

    for i in range(Lwq):
        temp = res[:, i]
        temp = list(np.squeeze(temp))
        temp = reduce(temp, q, centered=True)
        # print(temp)
        ans.append(temp)

    # Now each vector in ans is an element of R, and has as many coefficients as
    # we'd expect.

    return ans

# q = 8
# w = 2
# print(BitDecomp([-3, 1, -1]))

def pow_of_2(x):
    #global lwq
    global w
    lwq =  math.floor(math.log(q) / math.log(w)) + 2

    z = x
    y = []
    for i in range(lwq):
        z = [j*(w**i) for j in x]
        #print(i,z,reduce(z))
        z=reduce(z, q, centered=True)
        if z==[0.0]:
            z=z*len(x)
        y.append(z)
    return y


# print(pow_of_2([1,2]))

def LHE_ParamsGen(Lambda=987654321):
    """
    Given security parameter lambda, initialise d, q, t, ChiKey,
    ChiErr, and w, where w is an integer > 1
    """
    global d
    global q
    global t
    global B
    global n

    genProbabilities()

    # chi_key = ChiKey(n)
    # chi_err = ChiErr(n)

    d = 16
    w = sample(B)

    return d, q, t


def LHE_KeyGen():
    global d

    f_prime, g, f_inv, h, f = Basic_KeyGen(d)
    Lwq = int(np.floor(np.log(q)/np.log(w)))+2

    e = []
    s = []

    for i in range(Lwq):
        e.append(reduce(ChiErr(n), q, centered=True))

    for i in range(Lwq):
        s.append(reduce(ChiErr(n), q, centered=True))


    temp = pow_of_2(f)

    # print("size of ")
    # print("s, ", str(len(s)))
    # print("h, ", str(len(h)))
    # print("e, ", str(len(e)))

    # print(temp)

    const = []
    for i in range(len(s)):
        const.append(reduce(np.polymul(h, s[i]), q, centered=True))

    const1 = []
    for i in range(len(s)):
        const1.append(reduce(np.polyadd(e[i], const[i]), q, centered=True))

    #const = reduce(np.polymul(h, s), q, centered=True)
    # print("const, ", str(len(const)))

    Gamma = []
    for i in range(len(s)):
        Gamma.append(reduce(np.polyadd(temp[i], const1[i]), q, centered=True))
    # print(const)
    #for x in temp:
    #    Gamma.append(reduce(x+const, q, centered=True))

    # print(len(Gamma))
    # print(Lwq)

    return h, f, Gamma


# LHE_ParamsGen()
# print(LHE_KeyGen())


def LHE_Encrypt(h, msg):
    # Same as Basic Encrypt
    global q
    global t
    global n

    e = reduce(ChiErr(n), q, centered=True)
    s = reduce(ChiErr(n), q, centered=True)


    # print("e = ", e)
    # print("s = ", s)

    delta = math.floor(q/t)
    c = [delta*x for x in msg]
    c = np.polyadd(c, e)
    H = reduce(np.polymul(h, s), q, centered=True)
    c = np.polyadd(c, H)
    return reduce(c, q, centered=True)

def LHE_Decrypt(f, c):
    # Same as Basic Decrypt

    global q
    global t

    # print(len(c))
    # print(len(f))
    # print(len(c[0]))
    print('c = ',c)

    M = reduce(np.polymul(f, c), q, centered=True)
    print('M = ',M)

    M = [round(x*t/q) for x in M]
    print('M = ',M)

    return reduce(M, t, centered=True)

def LHE_Addition(c1, c2):

    global q

    A = np.polyadd(c1, c2)
    return reduce(A, q, centered=True)


def LHE_KeySwitch(cmult, evk):
    Dqw = BitDecomp(cmult)

    res = []
    Lwq = int(np.floor(np.log(q)/np.log(w)))+2

    for i in range(Lwq):
        temp = np.polymul(Dqw[i], evk[i])
        temp = reduce(temp, q, centered=True)
        res.append(temp)

    # temp = np.polymul(Dqw[i], evk[i])
    # res = reduce(temp, q, centered=True)

    return res


def LHE_Multiply(c1, c2, evk):
    global t
    global q

    const = t/q
    temp = np.polymul(c1, c2)

    temp = [np.round(x*const) for x in temp]

    cmultbar = reduce(temp, q, centered=True)
    temp = LHE_KeySwitch(cmultbar, evk)

    val = np.zeros(n)

    for x in temp:
        val = np.polyadd(val, x)

    val = reduce(val, q, centered=True)

    # print(val)

    return val
















LHE_ParamsGen()
h, f, Gamma = LHE_KeyGen()
# print('h = ',h)
# print('f = ',f)
# print('Gamma = ',Gamma)
# print('n = ',n)
# print('q = ',q)
# print('t = ',t)

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

c1 = LHE_Encrypt(h, reduce(msg1, t, centered=True))
c2 = LHE_Encrypt(h, reduce(msg2, t, centered=True))

final_msg = LHE_Decrypt(f, LHE_Addition(c1, c2))
print("Final msg is ", final_msg)

print("First cipher text is: ", c1)
print("Second cipher text is: ", c2)
print("Homomorphic Addition is ", reduce(final_msg, t))

final_msg = LHE_Decrypt(f, LHE_Multiply(c1, c2, Gamma))
print("Final msg is ", final_msg)

final_msg = reduce(final_msg, q, centered=True)

print("Final msg is ", final_msg)
print("Multiplication is ", reduce(final_msg, t))
