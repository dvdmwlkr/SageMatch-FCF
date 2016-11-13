︠90f3acbf-7a5c-4b45-a83e-58a582667ae7︠
sage_server.MAX_OUTPUT = 1000000

import time

vers = version()
print vers; print

def nearest(x):
    return floor(x+0.5)

def PrimeDivisorCount(n):
    pd = prime_divisors(n)
    return len(pd)

#F = (x for x in range(15,5001000,2) if (PrimeDivisorCount(x) > 1) and not is_square(x) and not is_prime(x) and (x %3 != 0) and (x%5 != 0))
F = (x for x in range(77,20050,2) if not is_square(x) and not is_prime(x) and (x %3 != 0) and (x%5 != 0))
#F = (x for x in range(343,345,2) if not is_square(x) and not is_prime(x) and (x %3 != 0) and (x%5 != 0))
#F = (x for x in range(100000000000000001,100000000000000100,2) if (PrimeDivisorCount(x) > 1) and not is_square(x) and not is_prime(x) and (x %3 != 0) and (x%5 != 0))
Failures = []
FailuresEvenLen = []
FailuresOddLen = []
MidDict = dict()

start_cpu = time.clock(); start_time = time.time()

xx = [1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000 ]
ix2 = 0
last_n = 15
countSoFar = 0 


found = 0
for n in F:
    #print "N:", n
    countSoFar += 1
    sq = floor(sqrt(n))
    A = [0, sq]
    B = [0, sq]
    C = [1, n - sq*sq]
    P = [1,sq]
    Q = [0, 1]
    i = 1
    #print"i: {0:4}  A: {1:4}  B: {2:4}  C: {3:4}  P: {4:4}  Q: {5:4}".format(i, A[i], B[i], C[i], P[i], Q[i])

    foundAtLeastOnefactor = False
    while C[i] != 1:
        k = i-1
        j = i
        i += 1
        A.append(floor((sq+B[j])/C[j]))
        B.append(A[i]*C[j]-B[j])
        C.append(C[k]+A[i]*(B[j]-B[i]))
        P.append((P[k]+A[i]*P[j]) % n)
        Q.append(Q[k]+A[i]*Q[j])
        #print"i: {0:4}  A: {1:4}  B: {2:4}  C: {3:4}  P: {4:4}  Q: {5:4}".format(i, A[i], B[i], C[i], P[i], Q[i])
        if (B[i] == B[i-1]) and C[i-1] > 2:
            if C[i-1] % 2 == 0:
                factr = C[i-1]/2
            else:
                factr = C[i-1]
            #print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  B[{0}]: {1}   factor: {2}".format(i, B[i], factr)
        foundFactor = False
        #if is_square(C[i]):
        if C[i] == 1:
            x = sqrt(C[i])
            x1 = [P[i]+x, P[i]-x]
            for x2 in x1:
                hcf = GCD(x2 ,n)
                if hcf > 1 and hcf < n:
                    foundFactor = True
                    foundAtLeastOnefactor = True
                    #print"{0} is a factor of {1}\n".format(hcf, n)
            if foundFactor:
                #print
                break;
        if C[i] == 1 and not foundAtLeastOnefactor:
            #print "Method failed - Sage factorisation: {0} = {1}\n".format(n, factor(n))
            #Failures.append(n)
            if len(A) % 2 == 1:
                FailuresEvenLen.append(n)
            else:
                FailuresOddLen.append(n)
            if len(A) % 2 == 1:
                ix = (len(A)+1)/2
                MidDict[n] = (B[ix], C[ix])
    if n > xx[ix2] and last_n < xx[ix2]:
        now_cpu = time.clock(); now_time = time.time()
        failures = len(FailuresEvenLen) + len(FailuresOddLen)
        print "N: {0:7d}   count: {1:6d}   failuresEvenLen: {2:5d}({3:.2f}%)   failuresOddLen: {4:5d}({5:.2f}%)   Success: {6:.2f}%   time: {7:.2f}".\
            format(n, countSoFar, len(FailuresEvenLen), float(len(FailuresEvenLen)*100/countSoFar), len(FailuresOddLen), float(len(FailuresOddLen)*100/countSoFar), float(((countSoFar-failures)*100)/countSoFar), now_time-start_time)
        ix2 += 1
    last_n = n
    if foundAtLeastOnefactor:
        found += 1
print "found:", found
︡d624b53e-f0f6-45ef-9616-bbefabca68de︡{"done":true,"error":"maximum time (=30000ms) exceeded - last error true"}︡









