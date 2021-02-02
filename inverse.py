import numpy as np

d = 4
phid = [1, 0, 1]
polynomials = []
q = 13

def reduce(p):
    global phid
    p[:] = [x % q for x in p]
    # print(p)
    Q, R = np.polydiv(p, phid)
    return R

def recurse(pos, d, q, myPol):
    global polynomials
    if pos == d+1:
        # print(myPol)
        newList = [int(x) for x in myPol]
        newList = reduce(newList)
        polynomials.append(newList)
        return
    for i in range (0, q):
        myPol.append(i)
        recurse(pos+1, d, q, myPol)
        myPol.pop()

def genPolynomials(q, d):
    myPol = []
    recurse(0, d, q, myPol)

genPolynomials(q, d/2)
print(len(polynomials))

print("Enter polynomial to find inverse")

poly = []
for i in range (0, 3):
    ele = int(input())
    poly.append(ele)

for p in polynomials:
    # print(p)
    f_f_inv = (reduce(np.polymul(p, poly))).tolist()
    if f_f_inv == [1]:
        print(p)
        break
