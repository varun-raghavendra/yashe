from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_gcdex
import random
import numpy as np
import math
import itertools

phid = np.zeros(4097)

phid[0]=1
phid[-1]=1
# print(phid)
# print(len(phid))
# phid = [1, 0, 0, 0, 0, 0, 0, 0, 0, 1]
d = 8192
# d = 16
q = math.pow(2, 126)+7
t = 256
w = math.pow(2, 32)
n = 4096
# n = 8
key_prob = []
err_prob = []
Berr = 48
Bkey = 1

sigma_err = Berr/6
sigma_key = Bkey/6

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

    # print("f_prime = ", f_prime)
    # print("g = ", g)
    # print("f = ", f)
    f_inv = []
    while True:
        try:
            X = (inverse(f))

        except:
            #print("Tried but no")
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
    global Berr
    global Bkey
    global n
    global key_prob

    genProbabilities()
    # print("Size of distributions\nChiKey ", len(key_prob))
    # print("ChiErr ", len(err_prob))
    # chi_key = ChiKey(n)
    # chi_err = ChiErr(n)

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

    # f = [-100.0, 0.0, 0.0, -50.0, 50.0, 0.0, 50.0, 1.0]
    # h = [264868.0, -88752.0, -245882.0, -405133.0, 216099.0, -500317.0, -422207.0, -19511.0]
    #
    # e = [[2.0, 1.0, 1.0, 1.0, 1.0, -1.0, 2.0, 1.0], [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, -1.0], [1.0, 0.0, 0.0, 0.0, 1.0, -2.0, 0.0, 1.0], [2.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, -1.0], [2.0, 1.0, -1.0, 1.0, 0.0, 1.0, -1.0, -2.0], [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, -4.0], [1.0, 1.0, 0.0, 1.0, 1.0, 0.0], [-2.0, -2.0, 2.0, -2.0, 1.0, 0.0, 1.0, -1.0], [2.0, 1.0, -1.0, 3.0, -1.0, 0.0, 1.0, -1.0], [-1.0, -1.0, 1.0, 1.0, 0.0, -1.0, 0.0, -1.0], [1.0, 1.0, -2.0, 0.0, -1.0, 1.0, 0.0, -1.0]]
    #
    # s = [[1.0, -1.0, -1.0, 1.0, 0.0, 1.0, -1.0, -2.0], [-2.0, 1.0, 0.0, -2.0, 0.0, 1.0, 0.0, 2.0], [2.0, -1.0, 0.0, -3.0], [-2.0, -2.0, -1.0, 1.0, -1.0, 1.0, 0.0, 2.0], [-1.0, 1.0, 2.0, 0.0, 0.0, 0.0, 2.0, 1.0], [-2.0, 0.0, 1.0, 1.0, -1.0, 1.0, -1.0, 3.0], [-2.0, -2.0, -1.0, 1.0, -3.0, 0.0, -2.0], [-1.0, 0.0, 1.0, 0.0, -1.0, -1.0, 1.0], [1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -2.0], [-3.0, 0.0, -1.0, 1.0, 0.0, -3.0, 0.0], [2.0, -2.0, 0.0, 0.0, 0.0, 0.0]]
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
    # print('f : ',f)
    # print('h : ',h)
    # print('e : ',e)
    # print('s : ',s)
    # print('Gamma : ',Gamma)

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








LHE_ParamsGen()
h, f, Gamma = LHE_KeyGen()


def demo():

    choice = None
    c1 = None
    c2 = None
    m1 = None
    m2 = None
    delta = np.floor(q/t)
    red = reduce([q], t, False)

    while(1):
        while(choice==None):
            print("\n\n****************************")
            print("1. Take input")
            print("2. Perform addition and print results")
            print("3. Perform multiplication")
            print("4. Print ciphertexts")
            print("5. Make large random inputs")
            print("6. Levels of multiplication")
            print("9. Exit")
            choice = int(input("Enter a choice: "))

        if(choice==1):
            print("Enter the plaintexts m1 and m2")
            s1 = input()
            s2 = input()
            m1 = s1.split()
            m2 = s2.split()
            m1 = [int(x) for x in m1]
            m2 = [int(x) for x in m2]
            print("The messages entered are:")
            print("message1: ", m1)
            print("message2: ", m2)
            c1 = LHE_Encrypt(h, reduce(m1, t, centered=True))
            c2 = LHE_Encrypt(h, reduce(m2, t, centered=True))
            choice=None
        elif(choice==2):
            if(m1==None or m2==None):
                continue

            final_msg = reduce(LHE_Decrypt(f, LHE_Addition(c1, c2)), t)
            expected_msg = reduce(np.polyadd(m1, m2), t)

            print("Homomorphic Addition is ", final_msg)
            print("Expected result is ", expected_msg)

            print("\n\n Are they the same?: ", bool(final_msg==expected_msg))

            print("Addition max allowed noise: ", (delta-red)/2)

            print("Addition 1 noise", norm(add_noise(f, c1, m1)))
            print("Addition 2 noise", norm(add_noise(f, c2, m2)))
            choice=None


        elif(choice==3):
            if(m1==None or m2==None):
                continue



            V = min(max(norm(add_noise(f, c1, m1)), norm(add_noise(f, c2, m2))), delta/2 - 1)

            Lwq = int(np.floor(np.log(q)/np.log(w)))+2
            const = n*t*(3+n*t*Bkey)*V+0.5*n*n*t*t*Bkey*(Bkey+t)+n*n*t*Lwq*w*Berr*Bkey

            c1c2, cmultbar = LHE_Multiply(c1, c2, Gamma)
            final_msg = LHE_Decrypt(f, c1c2)
            final_msg = reduce(final_msg, t, centered=True)
            expected_msg = reduce(np.polymul(m1, m2), t)
            print("Multiplication is ", final_msg)
            print("Expected : ", expected_msg)
            print("\n\n Are they the same?: ", bool(final_msg==expected_msg))
            mulnoise = norm(mul_noise(f, cmultbar, m1, m2))
            print("Multiplication noise", mulnoise)
            print("Multiplication max allowed noise: ", const)
            choice=None

        elif(choice==4):
            print("First cipher text is: ", c1)
            print("Second cipher text is: ", c2)

            choice=None

        elif(choice==5):
            length = np.random.randint(20)+1
            m1 = np.random.randint(30, size=(length))
            m2 = np.random.randint(30, size=(length))

            m1 = reduce(m1, t, True)
            m2 = reduce(m2, t, True)

            print("The messages generated are:")
            print("message1: ", m1)
            print("message2: ", m2)
            c1 = LHE_Encrypt(h, reduce(m1, t, centered=True))
            c2 = LHE_Encrypt(h, reduce(m2, t, centered=True))

            choice=None

        elif(choice==6):
            if(m1==None or m2==None):
                continue

            # print('Enter 3rd message :')
            # s3 = input()
            # m3 = s3.split()
            # m3 = [int(x) for x in m3]

            c1 = LHE_Encrypt(h, reduce(m1, t, centered=True))
            c2 = LHE_Encrypt(h, reduce(m2, t, centered=True))
            #c3 = LHE_Encrypt(h, reduce(m3, t, centered=True))

            V = min(max(norm(add_noise(f, c1, m1)), norm(add_noise(f, c2, m2))), delta/2 - 1)

            Lwq = int(np.floor(np.log(q)/np.log(w)))+2
            const = n*t*(3+n*t*Bkey)*V+0.5*n*n*t*t*Bkey*(Bkey+t)+n*n*t*Lwq*w*Berr*Bkey

            c1c2, cmultbar = LHE_Multiply(c1, c2, Gamma)
            final_msg = LHE_Decrypt(f, c1c2)
            final_msg = reduce(final_msg, t, centered=True)
            print("Multiplication in level 1 is ", reduce(final_msg, t))
            print("Expected : ", reduce(np.polymul(m1,m2),t))
            diff = np.polysub(reduce(final_msg, t), reduce(np.polymul(m1, m2), t))
            print('Difference : ', diff)
            comparison = diff==[0]
            #cc3, c3multbar = LHE_Multiply(c1c2, c3, Gamma)
            #new_final = LHE_Decrypt(f, cc3)
            #new_final = reduce(new_final, t, centered=True)
            #print('Multiplication level 3 is ', reduce(new_final, t))
            #print('Expected : ',reduce(np.polymul(m3, np.polymul(m1, m2)), t))

            mulnoise = norm(mul_noise(f, cmultbar, m1, m2))
            print("Multiplication noise", mulnoise)
            #print('Multiplication noise level 2 :', norm(mul_noise(f, c3multbar, m3, final_msg)))
            print("Multiplication max allowed noise: ", const)

            old_final = final_msg
            cc_prev = c1c2
            iterator = 2

            

            while(comparison.all()):
                try:
                    print('Enter next message :')
                    s3 = input()
                    m3 = s3.split()
                    m3 = [int(x) for x in m3]
                    c3 = LHE_Encrypt(h, reduce(m3, t, centered=True))
                    ccc, cimultbar = LHE_Multiply(cc_prev, c3, Gamma)
                    new_final = LHE_Decrypt(f, ccc)
                    new_final = reduce(new_final, t, centered=True)
                    print('Multiplication in level',iterator,  'is ', reduce(new_final, t))
                    print('Expected : ',reduce(np.polymul(m3, old_final), t))
                    diff = np.polysub(reduce(new_final, t), reduce(np.polymul(m3, old_final), t))
                    print('Difference : ', diff)
                    print('Multiplication noise in level',iterator,' :', norm(mul_noise(f, cimultbar, m3, old_final)))
                    print("Multiplication max allowed noise: ", const)

                    cc_prev=ccc
                    old_final=new_final
                    comparison = diff==[0]
                    iterator += 1

                    print("Do u wish to repeat :(y/n)")
                    string = input()
                    if string=='y':
                        continue
                    else:
                        break
                except:
                    print('Exception occured')
                    break


            choice=None




        else:
            break

    print("END OF DEMO")

demo()
