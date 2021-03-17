from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_gcdex
import random
import numpy as np
import math
import itertools

# Strategy
# The global variables on the outside stay as is, and we add a 
# self to each function, at the start, to show that it's a class 
# function. Then, we refer to each "global" variable with self.
# Then, we make some functions into private by putting __. 

class lfhe:
    phid = np.zeros(4097)

    phid[0]=1
    phid[-1]=1
    d = 8192
    q = math.pow(2, 126)+7
    t = 256
    w = math.pow(2, 32)
    n = 4096
    key_prob = []
    err_prob = []
    Berr = 48
    Bkey = 1

    sigma_err = Berr/6
    sigma_key = Bkey/6

    def reduce(p, mod, centered=True, divide=True):
        # Put the given polynomial in the range [-mod/2, mod/2) if centered,
        # and to [0,mod) if not centered.
        phid = self.phid

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
        global key_prob
        global err_prob
        global Berr
        global Bkey
        sum = 0
        for i in range (-1*Bkey, Bkey+1):
            x = np.e**(-(np.pi)*(i**2)/sigma_key**2)
            key_prob.append(x)
            sum += x
        key_prob = [x / sum for x in key_prob]

        sum = 0
        for i in range (-1*Berr, Berr+1):
            x = np.e**(-(np.pi)*(i**2)/sigma_err**2)
            err_prob.append(x)
            sum += x
        err_prob = [x / sum for x in err_prob]

    def sample(B, prob):
        # Returns a random sample from -B+1 to B, according to the
        # probability space
        return np.random.choice(np.arange(-1*B, B+1), p=prob)

    def ChiErr(n):
        # Makes a small coefficient error polynomial from the Chi function
        global Berr
        global err_prob
        poly = []
        for i in range (0, n):
            ele = sample(Berr, err_prob)
            poly.append(ele)
        return poly

    def ChiKey(n):
        # Makes a small coefficient key polynomial from the Chi function
        global Bkey
        poly = []
        for i in range (0, n):
            ele = sample(Bkey, key_prob)
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

        f_prime = ChiKey(d)
        g = ChiKey(d)
        f = [x*t for x in f_prime]
        f[-1] += 1
        f = reduce(f, q, centered=True)

        f_inv = []
        while True:
            try:
                X = (inverse(f))

            except:
                f_prime = ChiErr(d)
                f = [x*t for x in f_prime]
                f[-1] += 1
                f = reduce(f, q, centered=True)
                continue
            f_inv = [int(x) for x in X]
            if f_inv != [-1]:
                break

        h = reduce(np.polymul(g, f_inv), q, centered=True)
        h_prime = [x*t for x in h]
        h = reduce(h_prime, q, centered=True)

        return f_prime, g, f_inv, reduce(h, q, centered=True), reduce(f, q, centered=True)


    def split_to_bits(x, w, Lwq):
        temp = []

        temp = np.zeros(Lwq)

        if(x<0):
            x = pow(w, Lwq)+x

        idx = 0

        while(x):
            t = x%w
            temp[idx]=t
            idx+=1
            x = x//w
        return temp


    def BitDecomp(x):
        Lwq = int(np.floor(np.log(q)/np.log(w)))+2

        res = []

        for x_coeff in x:
            temp = split_to_bits(x_coeff, w, Lwq)
            res.append(temp)


        # print(res)
        res = np.asarray(res)

        ans = []

        for i in range(Lwq):
            temp = res[:, i]
            temp = list(np.squeeze(temp))
            temp = reduce(temp, q, centered=True)
            # print(temp)
            ans.append(temp)

        return ans


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


    def LHE_ParamsGen(Lambda=987654321):
        """
        Given security parameter lambda, initialise d, q, t, ChiKey,
        ChiErr, and w, where w is an integer > 1
        """
        global d
        global q
        global t
        global Berr
        global Bkey
        global n
        global key_prob

        genProbabilities()
        print("Size of distributions\nChiKey ", len(key_prob))
        print("ChiErr ", len(err_prob))

        w = sample(Bkey, key_prob)

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


        const = []
        for i in range(len(s)):
            const.append(reduce(np.polymul(h, s[i]), q, centered=True))

        const1 = []
        for i in range(len(s)):
            const1.append(reduce(np.polyadd(e[i], const[i]), q, centered=True))


        Gamma = []
        for i in range(len(s)):
            Gamma.append(reduce(np.polyadd(temp[i], const1[i]), q, centered=True))

        return h, f, Gamma



    def LHE_Encrypt(h, msg):
        # Same as Basic Encrypt
        global q
        global t
        global n

        e = reduce(ChiErr(n), q, centered=True)
        s = reduce(ChiErr(n), q, centered=True)

        # e =  [-2.0, -2.0, -1.0, 2.0, 0.0, 0.0, -1.0, -4.0]
        # s =  [-1.0, -2.0, -2.0, -1.0, 1.0, 0.0, 0.0, 0.0]

        #print("e = ", e)
        #print("s = ", s)

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
        #print('c = ',c)

        M = reduce(np.polymul(f, c), q, centered=True)
        #print('M = ',M)

        M = [round(x*t/q) for x in M]
        #print('M = ',M)

        return reduce(M, t, centered=True)

    def LHE_Addition(c1, c2):

        global q

        A = np.polyadd(c1, c2)
        return reduce(A, q, centered=True)


    def LHE_KeySwitch(cmult, evk):
        global n
        Dqw = BitDecomp(cmult)

        res = []
        Lwq = int(np.floor(np.log(q)/np.log(w)))+2

        for i in range(Lwq):
            temp = np.polymul(Dqw[i], evk[i])
            temp = reduce(temp, q, centered=True)
            res.append(temp)

        sum = np.zeros(n)
        for i in range(Lwq):
            sum = np.polyadd(sum,res[i])

        # temp = np.polymul(Dqw[i], evk[i])
        sum = reduce(sum, q, centered=True)

        return sum


    def LHE_Multiply(c1, c2, evk):
        global t
        global q

        const = t/q
        temp = np.polymul(c1, c2)

        temp = [np.round(x*const) for x in temp]

        cmultbar = reduce(temp, q, centered=True)
        #print('cmultbar : ',cmultbar)
        temp = LHE_KeySwitch(cmultbar, evk)

        #val = np.zeros(n)

        # for x in temp:
        #     val = np.polyadd(val, x)

        val = temp
        val = reduce(val, q, centered=True)

        # print(val)

        return val, cmultbar




    def add_noise(f, c, m):
        delta = np.floor(q/t)
        temp = np.polymul(f,c)

        del_m = [delta*x for x  in m]

        v = np.polysub(temp, del_m)

        return v

    def mul_noise(f, cmult, m1, m2):
        delta = np.floor(q/t)
        # temp = np.polymul(f,f)
        temp = np.polymul(f, cmult)
        m1m2 = reduce(np.polymul(m1, m2), t, True)


        del_m = [delta*x for x  in m1m2]

        noise = np.polysub(temp, del_m)
        return noise



    def norm(c):
        return np.max(np.abs(c))










