︠a25f513e-cd03-4168-82a8-47e95403519cs︠
import time

print version(); print

def PrimeDivisorCount(n):
    pd = prime_divisors(n)
    return len(pd)

#F = (x for x in range(803,2200,2) if (PrimeDivisorCount(x) > 1) and not is_square(x) and not is_prime(x) and (x %3 != 0) and (x%5 != 0))
F = (x for x in range(803,20050,2) if not is_square(x) and not is_prime(x) and (x %3 != 0) and (x%5 != 0))
Failures = []
FailuresEvenLen = []
FailuresOddLen = []
MidDict = dict()

start_cpu = time.clock(); start_time = time.time()

xx = [1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000 ]
ix2 = 0
last_n = 15
countSoFar = 0 

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
        P.append(P[k]+A[i]*P[j])
        Q.append(Q[k]+A[i]*Q[j])
        #print"i: {0:4}  A: {1:4}  B: {2:4}  C: {3:4}  P: {4:4}  Q: {5:4}".format(i, A[i], B[i], C[i], P[i], Q[i])
        #if (B[i] == B[i-1]) and C[i-1] > 2:
        #    if C[i-1] % 2 == 0:
        #        factr = C[i-1]/2
        #    else:
        #        factr = C[i-1]
        #    print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  B[{0}]: {1}   factor: {2}".format(i-1, B[i], factr)
        foundFactor = False
        if is_square(C[i]):
            x = sqrt(C[i])
            x1 = [P[i]+x, P[i]-x]
            for x2 in x1:
                hcf = GCD(x2 ,n)
                if hcf > 1 and hcf < n:
                    foundFactor = True
                    foundAtLeastOnefactor = True
                    #print"{0} is a factor of {1}".format(hcf, n)
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
        if ix2 ==4:
            print "Even Failures({0}): \n  {1}   \nOdd Failures({2}): \n  {3}".format(len(FailuresEvenLen), FailuresEvenLen, len(FailuresOddLen), FailuresOddLen)
        ix2 += 1
    last_n = n
failuresPlus2 =[]
failuresMinus2 =[]
uvPlus2= []
uvMinus2= []
for k in FailuresEvenLen:
     s = QuadraticField(k).gen()
     cf = continued_fraction(s)
     #print '%2d %d %s'%(k, len(cf.period()), cf)
     for x in range(len(cf.period())+1):
         if cf.p(x)^2 - k*cf.q(x)^2 == 1:
             halfWay = cf.p((x-1)/2)^2 - k*cf.q((x-1)/2)^2
             #print "n:", k, halfWay
             if halfWay == 2:
                 failuresPlus2.append(k)
                 uvPlus2.append((cf.p((x-1)/2), cf.q((x-1)/2)))
             else:
                 failuresMinus2.append(k)
                 uvMinus2.append((cf.p((x-1)/2), cf.q((x-1)/2)))
             break
         #print x, cf.p(x), cf.q(x), cf.p(x)^2 - k*cf.q(x)^2
#print "failuresPlus2({0}):\n {1}".format(len(failuresPlus2),failuresPlus2)
#print "uvPlus2:\n{0}".format(uvPlus2)
#print "failuresMinus2({0}):\n {1}".format(len(failuresMinus2),failuresMinus2)
#print "uvMinus2:\n{0}\n".format(uvMinus2)

mid_cpu = time.clock(); mid_time = time.time()
print "Time to split = {0:.2f}".format(mid_time-start_time)

class qfm2:
    """Quadratic Field Class (-2_)"""

    def __init__(self, part1, part2):
        self.p1 = part1
        self.p2 = part2

    def __add__(self, other):
        return qfm2(self.p1+other.p1, self.p2+other.p2)

    def __sub__(self, other):
        return qfm2(self.p1-other.p1, self.p2-other.p2)

    def __mul__(self, other):
        return qfm2(self.p1*other.p1-2*self.p2*other.p2, self.p2*other.p1+self.p1*other.p2)

    def __div__(self, other):
        q = other.p1^2+2*other.p2^2
        return qfm2((self.p1*other.p1+2*self.p2*other.p2)/q, (self.p2*other.p1-self.p1*other.p2)/q)

    def __neg__(self):
        return qfm2(-self.p1, -self.p2)

    def __eq__(self, other):
        return (self.p1 == other.p1 and self.p2 == other.p2)
    def __ne__(self, other):
        return (self.p1 != other.p1 or self.p2 != other.p2)

    def __repr__(self):
        return "({0}, {1})".format(self.p1, self.p2)

    def sqrt(self):
        x = sqrt((self.p1 + sqrt(self.p1*self.p1 + 2*self.p2*self.p2))/2)
        y = sqrt(self.p1/2) if x == 0 else self.p2/(2*x)
        return qfm2(x, y)

    def ni(self):
        return qfm2(floor(self.p1+0.5), floor(self.p2+0.5))

    def __mod__(self, other):
        return self - other*((self/other).ni())

def norm(self):
    return self.p1*self.p1+2*self.p2*self.p2

def is_qfm2square(self):
    def is_int(x):
        return (x - floor(x) == 0)
    c = sqrt(self.p1^2+2*self.p2^2)
    x = sqrt((self.p1+c)/2)
    y = sqrt(-self.p1/2) if x == 0 else self.p2/(2*x)
    if is_int(y) and ((is_int(c) and is_int(x) and x!= 0) or (x == 0)):
        return True
    else:
        return False

def GCD_qfm2(self, other):
    j,k = self, other
    if norm(j) < norm(k):
        j,k = k,j
    #print "j", j, "k:", k
    while true:
        q = j/k
        qq = q.ni()
        #print "j/k", q
        remainder = j - k*qq
        #print "qq", qq, "remainder", remainder, "k:", k
        if remainder == qfm2(0,0):
            break
        j,k = k,remainder
        #print "j", j, "k:", k
    return qfm2(abs(k.p1), abs(k.p2))

fl = []
success_count = failure_count = 0
for m in enumerate(failuresMinus2):
    #print "{0}    {1}".format(m[1], uvMinus2[m[0]][0])
    a = qfm2(m[1], 0)
    b = qfm2(uvMinus2[m[0]][0], 1)
    c = GCD_qfm2(a,b)
    print "\nN:", m[1]
    #print "a,b,c",a,b,c

    x1 = float(sqrt((c.p1 + float(sqrt(m[1])))/2))
    x2 = float(abs(c.p2/(2*x1)))
    sq = qfm2(x1,x2)
    rounded = sq.ni()
    #print "{0}^2 - {1}*{2}^2 = -2".format(b.p1, m[1], uvMinus2[m[0]][1])
    #print "HCF(({0},0),({1},{2})) = ({3},{4})      sq: {5}  rounded: {6} rounded^2: {7}\n".format(a, b.p1, b.p2, c.p1, c.p2, sq, rounded, rounded*rounded)
    A = [0, rounded]
    B = [0, rounded]
    C = [qfm2(1,0), c - rounded*rounded]
    P = [qfm2(1,0), rounded]
    Q = [qfm2(0,0), qfm2(1,0)]
    i = 1
    #print "A: {0:<8}   B: {1:<8}   C: {2:<8}   P: {3:<10}  i: {4}".format(A[-1], B[-1], C[-1], P[-1], i)
    for zz in range(200):
        k = i-1
        j = i
        i += 1
        tmp = (sq+B[j])/C[j]
        A.append(((sq+B[j])/C[j]).ni())
        B.append(A[i]*C[j]-B[j])
        C.append(C[k]+A[i]*(B[j]-B[i]))
        #P.append((P[k]+A[i]*P[j]))
        P.append((P[k]+A[i]*P[j]) % c)
        Q.append(Q[k]+A[i]*Q[j])
        #print "A: {0:<8}   B: {1:<8}   C: {2:<8}   P: {3:<10}  i: {4}".format(A[-1], B[-1], C[-1], P[-1], i)
        ct = -C[-1] if i%2 > 0 else C[-1] 
        if is_qfm2square(ct):
            zero = qfm2(0,0)
            cs = sqrt(ct)
            str = "i: {0}  sqrt{1} = {2}".format(i, ct, cs)
            print str
            factor = []
            t1 = P[-1] + cs
            if t1 != zero:
                f1 = GCD_qfm2(t1, c)
                factor.append(norm(f1))
            t2 = P[-1] - cs
            if t2 != zero:
                f2 = GCD_qfm2(t2, c)
                factor.append(norm(f2))
#            if len(factor) == 2:
#                print "{0} - factors {1}".format(a.p1, factor)
            if len(factor) == 2 and factor[0]*factor[1] == a.p1:
                print "{0} = {1}*{2}".format(a.p1, factor[0], factor[1])
                success_count += 1
                break
            if C[-1] == qfm2(1,0):
                fl.append((m[1], str))
                failure_count += 1
                print "\n"
                break
print "\nsuccess_count: {0}    failure_count: {1}   failure percentage: {2:.2f}%\n".format(success_count, failure_count, float((failure_count * 100)/(success_count+failure_count)))
print "Failure list:"
for s in fl:
    print s

now_cpu = time.clock(); now_time = time.time()
print "Time to end = {0:.2f}".format(now_time-mid_time)
︡655fc138-2922-4ecd-b3ef-4ca35d86165c︡{"stdout":"SageMath version 7.3, Release Date: 2016-08-04\n\n"}︡{"stdout":"N:    1001   count:     22   failuresEvenLen:     2(9.09%)   failuresOddLen:     0(0.00%)   Success: 90.91%   time: 0.02\nN:    2009   count:    152   failuresEvenLen:    11(7.24%)   failuresOddLen:     2(1.32%)   Success: 91.45%   time: 0.08\nN:    5017   count:    580   failuresEvenLen:    40(6.90%)   failuresOddLen:    13(2.24%)   Success: 90.86%   time: 0.31"}︡{"stdout":"\nN:   10001   count:   1344   failuresEvenLen:    83(6.18%)   failuresOddLen:    26(1.93%)   Success: 91.89%   time: 0.73"}︡{"stdout":"\nN:   20003   count:   2967   failuresEvenLen:   166(5.59%)   failuresOddLen:    68(2.29%)   Success: 92.11%   time: 1.81"}︡{"stdout":"\nEven Failures(166): \n  [803, 959, 1003, 1127, 1331, 1343, 1411, 1519, 1679, 1687, 1691, 2023, 2123, 2147, 2263, 2599, 2603, 2651, 2747, 2807, 2863, 3239, 3247, 3403, 3479, 3587, 3647, 3667, 3707, 3791, 3827, 3971, 4151, 4223, 4247, 4319, 4487, 4499, 4579, 4763, 5027, 5327, 5339, 5371, 5383, 5543, 5627, 6167, 6191, 6239, 6319, 6347, 6443, 6499, 6503, 6559, 6611, 6727, 6859, 6887, 7223, 7339, 7343, 7403, 7567, 7571, 7619, 8023, 8227, 8279, 8459, 8483, 8519, 8531, 8651, 8687, 8899, 8927, 9023, 9079, 9247, 9407, 9799, 10199, 10291, 10307, 10327, 10423, 10447, 10727, 10747, 10823, 10951, 10999, 11023, 11099, 11207, 11263, 11303, 11327, 11543, 11711, 11723, 11879, 12167, 12179, 12319, 12431, 12587, 12767, 13067, 13123, 13223, 13379, 13391, 13439, 13739, 13823, 13943, 14003, 14191, 14267, 14507, 14603, 14611, 14903, 14959, 15127, 15479, 15623, 15967, 16019, 16199, 16379, 16439, 16571, 16639, 16643, 16847, 17063, 17083, 17111, 17311, 17363, 17639, 17687, 17867, 18071, 18227, 18599, 18631, 18767, 18779, 18983, 18991, 19103, 19307, 19439, 19547, 19567, 19723, 19771, 19871, 19879, 19883, 19967]   \nOdd Failures(68): \n  [1313, 1853, 2197, 2249, 2813, 3029, 3281, 3293, 3341, 3653, 3973, 4469, 4913, 5213, 5353, 5629, 5713, 5837, 5933, 6253, 6893, 8321, 8653, 9197, 9773, 9953, 10397, 11029, 11141, 11401, 11453, 11509, 11629, 11773, 11849, 12053, 12389, 12461, 12557, 12629, 12773, 13549, 13801, 14089, 14453, 15133, 15509, 15529, 15613, 15977, 16109, 16133, 16601, 16609, 16733, 17177, 17261, 17693, 17797, 17849, 17893, 18317, 18409, 18989, 19109, 19409, 19633, 19637]\n"}︡{"stdout":"Time to split = 2.83\n"}︡{"stdout":"\nN: 803\ni: 2  sqrt(-1, -2) = (1, -1)\n803 = 11*73\n\nN: 1003\ni: 6  sqrt(1, 0) = (1, 0)\n1003 = 59*17\n\nN: 1331\ni: 13  sqrt(-1, -2) = (1, -1)"}︡{"stdout":"\ni: 25  sqrt(-1, -2) = (1, -1)\ni: 38  sqrt(1, 0) = (1, 0)"}︡{"stdout":"\n\n\n\nN: 1411\ni: 14  sqrt(1, 0) = (1, 0)"}︡{"stdout":"\n1411 = 17*83\n\nN: 1691\ni: 6  sqrt(-1, -2) = (1, -1)\n1691 = 89*19\n\nN: 2123\ni: 8  sqrt(-1, -2) = (1, -1)\n2123 = 11*193\n\nN: 2147\ni: 7  sqrt(-7, 4) = (1, 2)"}︡{"stdout":"\ni: 11  sqrt(-1, -2) = (1, -1)\n2147 = 19*113\n\nN: 2603\ni: 3  sqrt(-1, -2) = (1, -1)\ni: 23  sqrt(-1, -2) = (1, -1)"}︡{"stdout":"\n2603 = 137*19\n\nN: 2651\ni: 12  sqrt(1, 0) = (1, 0)\n2651 = 241*11\n\nN: 2747\ni: 7  sqrt(-1, -2) = (1, -1)"}︡{"stdout":"\n2747 = 41*67\n\nN: 3403\ni: 14  sqrt(-1, 2) = (1, 1)"}︡{"stdout":"\ni: 26  sqrt(9, 0) = (3, 0)"}︡{"stdout":"\ni: 30  sqrt(-7, -4) = (1, -2)\ni: 32  sqrt(-1, -2) = (1, -1)\n3403 = 41*83\n\nN: 3587\ni: 2  sqrt(-1, -2) = (1, -1)\n3587 = 17*211\n\nN: 3667\ni: 2  sqrt(1, 0) = (1, 0)\n3667 = 193*19\n\nN: 3707\ni: 10  sqrt(-1, -2) = (1, -1)"}︡{"stdout":"\n3707 = 11*337\n\nN: 3827\ni: 5  sqrt(-7, -4) = (1, -2)\ni: 8  sqrt(-1, 2) = (1, 1)\ni: 14  sqrt(-1, 2) = (1, 1)\ni: 17  sqrt(-7, -4) = (1, -2)"}︡{"stdout":"\ni: 22  sqrt(1, 0) = (1, 0)\n\n\n\nN: 3971\ni: 3  sqrt(-1, -2) = (1, -1)\ni: 31  sqrt(-1, -2) = (1, -1)"}︡{"stdout":"\n3971 = 361*11\n\nN: 4499\ni: 6  sqrt(-1, -2) = (1, -1)\ni: 16  sqrt(-7, 4) = (1, 2)"}︡{"stdout":"\ni: 28  sqrt(-7, 4) = (1, 2)\n4499 = 11*409\n\nN: 4579\ni: 5  sqrt(7, -6) = (3, -1)"}︡{"stdout":"\ni: 20  sqrt(1, 0) = (1, 0)"}︡{"stdout":"\n4579 = 241*19\n\nN: 4763\ni: 10  sqrt(-1, -2) = (1, -1)\n4763 = 433*11\n\nN: 5027\ni: 8  sqrt(-1, -2) = (1, -1)"}︡{"stdout":"\n5027 = 11*457\n\nN: 5339\ni: 11  sqrt(-7, -4) = (1, -2)\ni: 13  sqrt(-1, 2) = (1, 1)\n5339 = 19*281\n\nN: 5371\ni: 4  sqrt(9, 0) = (3, 0)"}︡{"stdout":"\n5371 = 41*131\n\nN: 5627\ni: 13  sqrt(7, -6) = (3, -1)\n5627 = 17*331\n\nN: 6347\ni: 4  sqrt(1, 0) = (1, 0)"}︡{"stdout":"\n6347 = 11*577\n\nN: 6443\ni: 3  sqrt(-1, 2) = (1, 1)\ni: 9  sqrt(-7, -4) = (1, -2)\ni: 33  sqrt(-7, -4) = (1, -2)"}︡{"stdout":"\n6443 = 379*17\n\nN: 6499\ni: 6  sqrt(1, 0) = (1, 0)\n6499 = 97*67\n\nN: 6611\ni: 5  sqrt(-1, 2) = (1, 1)\n6611 = 601*11\n\nN: 6859\ni: 50  sqrt(1, 0) = (1, 0)"}︡{"stdout":"\n\n\n\nN: 7339\ni: 22  sqrt(1, 0) = (1, 0)"}︡{"stdout":"\n7339 = 41*179\n\nN: 7403\ni: 16  sqrt(-7, 4) = (1, 2)"}︡{"stdout":"\n7403 = 673*11\n\nN: 7571\ni: 7  sqrt(-7, 4) = (1, 2)\ni: 20  sqrt(-1, -2) = (1, -1)"}︡{"stdout":"\ni: 24  sqrt(-1, -2) = (1, -1)\ni: 37  sqrt(-7, 4) = (1, 2)"}︡{"stdout":"\ni: 44  sqrt(1, 0) = (1, 0)\n\n\n\nN: 7619\ni: 16  sqrt(1, 0) = (1, 0)"}︡{"stdout":"\n7619 = 19*401\n\nN: 8227\ni: 4  sqrt(9, 0) = (3, 0)\ni: 12  sqrt(-1, 2) = (1, 1)\n8227 = 19*433\n\nN: 8459\ni: 4  sqrt(-1, -2) = (1, -1)\n8459 = 769*11\n\nN: 8483\ni: 7  sqrt(-1, 2) = (1, 1)"}︡{"stdout":"\ni: 14  sqrt(-7, -4) = (1, -2)\ni: 33  sqrt(-1, 2) = (1, 1)"}︡{"stdout":"\n8483 = 499*17\n\nN: 8531\ni: 7  sqrt(-1, -2) = (1, -1)\n8531 = 449*19\n\nN: 8651\ni: 12  sqrt(-1, -2) = (1, -1)"}︡{"stdout":"\n8651 = 41*211\n\nN: 8899\ni: 3  sqrt(-1, 2) = (1, 1)\ni: 9  sqrt(-7, -4) = (1, -2)\ni: 21  sqrt(9, 0) = (3, 0)"}︡{"stdout":"\n8899 = 11*809\n\nN: 10291\ni: 12  sqrt(-1, 2) = (1, 1)\ni: 25  sqrt(-7, -4) = (1, -2)"}︡{"stdout":"\ni: 28  sqrt(-1, -2) = (1, -1)\n10291 = 251*41\n\nN: 10307\ni: 5  sqrt(-7, 4) = (1, 2)\ni: 10  sqrt(7, 6) = (3, 1)\ni: 28  sqrt(-1, -2) = (1, -1)"}︡{"stdout":"\ni: 34  sqrt(-1, -2) = (1, -1)\ni: 52  sqrt(7, 6) = (3, 1)"}︡{"stdout":"\ni: 57  sqrt(-7, 4) = (1, 2)\ni: 62  sqrt(1, 0) = (1, 0)\n\n\n\nN: 10747\ni: 6  sqrt(1, 12) = (3, 2)"}︡{"stdout":"\ni: 18  sqrt(7, 6) = (3, 1)\n10747 = 977*11\n\nN: 11099\ni: 8  sqrt(-1, 2) = (1, 1)"}︡{"stdout":"\n11099 = 11*1009\n\nN: 11723\ni: 11  sqrt(-7, -4) = (1, -2)\ni: 28  sqrt(1, 0) = (1, 0)"}︡{"stdout":"\n11723 = 617*19\n\nN: 12179\ni: 6  sqrt(1, 0) = (1, 0)\n12179 = 19*641\n\nN: 12587\ni: 4  sqrt(-7, -4) = (1, -2)\n12587 = 41*307\n\nN: 13067\ni: 12  sqrt(-1, -2) = (1, -1)"}︡{"stdout":"\n13067 = 179*73\n\nN: 13123\ni: 10  sqrt(1, 0) = (1, 0)\n13123 = 11*1193\n\nN: 13379\ni: 6  sqrt(1, 0) = (1, 0)"}︡{"stdout":"\n13379 = 787*17\n\nN: 13739\ni: 3  sqrt(-1, 2) = (1, 1)\n13739 = 11*1249\n\nN: 14003\ni: 8  sqrt(-1, -2) = (1, -1)\n14003 = 1273*11\n\nN: 14267\ni: 8  sqrt(-7, -4) = (1, -2)"}︡{"stdout":"\n14267 = 1297*11\n\nN: 14507\ni: 7  sqrt(-1, -2) = (1, -1)\n14507 = 89*163\n\nN: 14603\ni: 13  sqrt(-1, -2) = (1, -1)"}︡{"stdout":"\n14603 = 17*859\n\nN: 14611\ni: 24  sqrt(1, 0) = (1, 0)"}︡{"stdout":"\n14611 = 769*19\n\nN: 16019\ni: 22  sqrt(1, 0) = (1, 0)"}︡{"stdout":"\n16019 = 193*83\n\nN: 16379\ni: 12  sqrt(-7, 4) = (1, 2)"}︡{"stdout":"\n16379 = 11*1489\n\nN: 16571\ni: 3  sqrt(-7, -4) = (1, -2)\ni: 5  sqrt(-7, -4) = (1, -2)\n16571 = 227*73\n\nN: 16643\ni: 3  sqrt(-1, -2) = (1, -1)\n16643 = 11*1513\n\nN: 17083\ni: 11  sqrt(-7, -4) = (1, -2)"}︡{"stdout":"\ni: 19  sqrt(-1, -2) = (1, -1)\n17083 = 1553*11\n\nN: 17363\ni: 8  sqrt(7, -6) = (3, -1)"}︡{"stdout":"\ni: 16  sqrt(7, 6) = (3, 1)\ni: 21  sqrt(-1, 2) = (1, 1)"}︡{"stdout":"\ni: 32  sqrt(-7, -4) = (1, -2)"}︡{"stdout":"\n17363 = 97*179\n\nN: 17867\ni: 11  sqrt(-1, 2) = (1, 1)\ni: 21  sqrt(-7, -4) = (1, -2)"}︡{"stdout":"\ni: 26  sqrt(7, -6) = (3, -1)\ni: 32  sqrt(7, 6) = (3, 1)\ni: 62  sqrt(7, 6) = (3, 1)"}︡{"stdout":"\n17867 = 1051*17\n\nN: 18227\ni: 3  sqrt(-7, 4) = (1, 2)\n18227 = 11*1657\n\nN: 18779\ni: 4  sqrt(-7, 4) = (1, 2)\ni: 15  sqrt(7, -6) = (3, -1)"}︡{"stdout":"\n18779 = 89*211\n\nN: 19307\ni: 5  sqrt(-1, 2) = (1, 1)\ni: 12  sqrt(-7, -4) = (1, -2)\ni: 33  sqrt(7, -6) = (3, -1)"}︡{"stdout":"\ni: 37  sqrt(7, -6) = (3, -1)\n19307 = 43*449\n\nN: 19547\ni: 13  sqrt(-1, -2) = (1, -1)"}︡{"stdout":"\ni: 23  sqrt(1, 12) = (3, 2)"}︡{"stdout":"\n19547 = 1777*11\n\nN: 19723\ni: 4  sqrt(7, 6) = (3, 1)\ni: 8  sqrt(7, -6) = (3, -1)\n19723 = 1793*11\n\nN: 19771\ni: 5  sqrt(9, 0) = (3, 0)"}︡{"stdout":"\n19771 = 1163*17\n\nN: 19883\ni: 4  sqrt(7, -6) = (3, -1)\ni: 16  sqrt(-1, -2) = (1, -1)"}︡{"stdout":"\n19883 = 337*59\n"}︡{"stdout":"\nsuccess_count: 63    failure_count: 5   failure percentage: 7.35%\n\n"}︡{"stdout":"Failure list:\n"}︡{"stdout":"(1331, 'i: 38  sqrt(1, 0) = (1, 0)')\n(3827, 'i: 22  sqrt(1, 0) = (1, 0)')\n(6859, 'i: 50  sqrt(1, 0) = (1, 0)')\n(7571, 'i: 44  sqrt(1, 0) = (1, 0)')\n(10307, 'i: 62  sqrt(1, 0) = (1, 0)')\n"}︡{"stdout":"Time to end = 8.63\n"}︡{"done":true}︡









