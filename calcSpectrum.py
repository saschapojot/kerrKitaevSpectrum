from consts import *
from multiprocessing import Pool
import matplotlib.pyplot as plt
# values of Hamiltonian matrix elements




def h2(k):
    '''

    :param k: 0,1,...,N-1
    :return:
    '''
    rst = -d1 * np.sin(dk * k)
    return rst


def h3(k, deltaK):
    rst = lmd / 2 * deltaK - mu1 / 2 - t1 * np.cos(dk * k)
    return rst



def coef3(k):
    rst = -lmd * (mu1 / 2 + t1 * np.cos(dk * k))
    return rst


def coef2(k):
    rst = d1 ** 2 * (np.sin(k * dk)) ** 2 + (mu1 / 2 + t1 * np.cos(k * dk)) ** 2 - lmd ** 2 * rho ** 2 / 4
    return rst


def coef1(k):
    rst = lmd * rho ** 2 * (mu1 / 2 + t1 * np.cos(k * dk))
    return rst


def coef0(k):
    rst = -rho ** 2 * (mu1 / 2 + t1 * np.cos(dk * k)) ** 2
    return rst


def solve2(k):
    '''
    for lambda=0
    :param k:
    :return:
    '''
    coef2Val = coef2(k)
    coef1Val = coef1(k)
    coef0Val = coef0(k)
    inCoefs = [coef2Val, coef1Val, coef0Val]
    return np.roots(inCoefs)


def solve4(k):
    '''
    for lambda != 0
    :param k:
    :return:
    '''

    coef3Val = coef3(k)
    coef2Val = coef2(k)
    coef1Val = coef1(k)
    coef0Val = coef0(k)
    inCoefs = [coef4Val, coef3Val, coef2Val, coef1Val, coef0Val]
    return np.roots(inCoefs)


def solveEqn(k):
    if np.abs(lmd) > tol:
        return [k, solve4(k)]
    else:
        return [k, solve2(k)]


pool1 = Pool(threadNum)
kAndRoot = pool1.map(solveEqn, kIndHalf)
pool1.close()
pool1.join()
for item in kAndRoot:
    kVal = item[0]
    rts = item[1]
    tmp = []
    for dkVal in rts:
        if np.abs(np.imag(dkVal)) < tol:
            tmp.append(np.real(dkVal))
    if len(tmp) > 0:
        deltaKAll[kVal] = tmp

for kVal in deltaKAll:
    dkVals=deltaKAll[kVal]
    itemTmp=[]
    for oneDeltakVal in dkVals:
        h3Val=h3(kVal,oneDeltakVal)
        Etmp=h0Val+h3Val*rho/oneDeltakVal
        itemTmp.append(Etmp)
    EkAll[kVal]=itemTmp


def h2L(k):
    rst=-d1*np.sin(dk*k)
    return rst
def h3L(k):
    rst=-mu1/2-t1*np.cos(k*dk)
    return rst
def EL(k):
    h2LVal=h2L(k)
    h3LVal=h3L(k)
    tmp=np.sqrt(h2LVal**2+h3LVal**2)
    return [k,h0LVal-tmp,h0LVal+tmp]

pool2=Pool(threadNum)
kE2ValsAll=pool2.map(EL,kIndHalf)
pool2.close()
pool2.join()


plt.figure()
#nonlin spectrum
for kVal in EkAll:
    EValsTmp=EkAll[kVal]
    kTrue=kVal*dk
    for E in EValsTmp:
        blk=plt.scatter(kTrue,E,c="black",s=1)

#lin spectrum
for item in kE2ValsAll:
    kVal=item[0]
    Em=item[1]
    Ep=item[2]
    kTrue=kVal*dk
    rd=plt.scatter(kTrue,Em,c="red",s=1)
    rd=plt.scatter(kTrue,Ep,c="red",s=1)
plt.legend([rd,blk],["linear","nonlinear"])
plt.xlabel("momentum")
plt.ylabel("energy")
plt.title("$\mu_{1}=$"+str(mu1)+", $t_{1}=$"+str(t1)+", $\Delta_{1}=$"+str(d1)+", $\lambda=$"+str(lmd))


outFileName="./spectrum2/specmu1"+str(mu1)+"t1"+str(t1)+"d1"+str(d1)+"lambda"+str(lmd)+".png"
plt.savefig(outFileName)
plt.close()



