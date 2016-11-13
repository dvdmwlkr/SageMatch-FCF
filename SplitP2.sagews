︠bd4f7061-afc5-417c-be90-fe02c67872d7︠
import time
from abc import ABCMeta, abstractmethod

print version(); print

def PrimeDivisorCount(n):
    pd = prime_divisors(n)
    return len(pd)

F = (x for x in range(119,20050,2) if not is_square(x) and not is_prime(x) and (x %3 != 0) and (x%5 != 0))
Failures = []
FailuresEvenLen = []
FailuresOddLen = []

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
             x1 = (x-1)/2
             halfWay = cf.p(x1)^2 - k*cf.q(x1)^2
             #print "n:", k, halfWay
             if halfWay == 2:
                 failuresPlus2.append(k)
                 uvPlus2.append((cf.p(x1), cf.q(x1)))
             else:
                 failuresMinus2.append(k)
                 uvMinus2.append((cf.p(x1), cf.q(x1)))
             break
         #print x, cf.p(x), cf.q(x), cf.p(x)^2 - k*cf.q(x)^2
print "failuresPlus2({0}):\n {1}".format(len(failuresPlus2),failuresPlus2)
print "uvPlus2:\n{0}".format(uvPlus2)
#print "failuresMinus2({0}):\n {1}".format(len(failuresMinus2),failuresMinus2)
#print "uvMinus2:\n{0}\n".format(uvMinus2)

mid_cpu = time.clock(); mid_time = time.time()
print "Time to split = {0:.2f}".format(mid_time-start_time)

class qf:
    """Quadratic Field Class - an abstract base class"""
    def __init__(self, part1, part2=0):
        self.p1 = part1
        self.p2 = part2
    def __eq__(self, other):
        return (self.p1 == other.p1 and self.p2 == other.p2)
    def __ne__(self, other):
        return (self.p1 != other.p1 or self.p2 != other.p2)
    def __repr__(self):
        return "({0}, {1})".format(self.p1, self.p2)
    def __mod__(self, other):
        return self - other*((self/other).ni())
    def qftuple(self):
        return self.p1,self.p2
    # The following methods must be implemented in any subclass
    __metaclass__=ABCMeta
    @abstractmethod
    def ni(self):
        pass
    @abstractmethod
    def __neg__(self):
        pass
    @abstractmethod
    def __sub__(self, other):
        pass
    @abstractmethod
    def __add__(self, other):
        pass
    @abstractmethod
    def __mul__(self, other):
        pass
    @abstractmethod
    def __div__(self, other):
        pass
    @abstractmethod
    def sqrt(self):
        pass
    @abstractmethod
    def norm(self):
        pass

class qfm2(qf):
    """Quadratic Field Class (-2)"""
    def ni(self):
        return qfm2(floor(self.p1+0.5), floor(self.p2+0.5))
    def __neg__(self):
        return qfm2(-self.p1, -self.p2)
    def __sub__(self, other):
        return qfm2(self.p1-other.p1, self.p2-other.p2)
    def __add__(self, other):
        return qfm2(self.p1+other.p1, self.p2+other.p2)
    def __mul__(self, other):
        return qfm2(self.p1*other.p1-2*self.p2*other.p2, self.p2*other.p1+self.p1*other.p2)
    def __div__(self, other):
        q = other.p1^2+2*other.p2^2
        return qfm2((self.p1*other.p1+2*self.p2*other.p2)/q, (self.p2*other.p1-self.p1*other.p2)/q)
    def sqrt(self):
        x = sqrt((self.p1 + sqrt(self.p1*self.p1 + 2*self.p2*self.p2))/2)
        y = sqrt(self.p1/2) if x == 0 else self.p2/(2*x)
        return qfm2(x, y)
    def norm(self):
        return self.p1*self.p1+2*self.p2*self.p2

class qfp2(qf):
    """Quadratic Field Class (+2)"""
    def ni(self):
        return qfp2(floor(self.p1+0.5), floor(self.p2+0.5))
    def __neg__(self):
        return qfp2(-self.p1, -self.p2)
    def __sub__(self, other):
        return qfp2(self.p1-other.p1, self.p2-other.p2)
    def __add__(self, other):
        return qfp2(self.p1+other.p1, self.p2+other.p2)
    def __mul__(self, other):
        return qfp2(self.p1*other.p1-2*self.p2*other.p2, self.p2*other.p1+self.p1*other.p2)
    def __mul__(self, other):
        return qfp2(self.p1*other.p1+2*self.p2*other.p2, self.p2*other.p1+self.p1*other.p2)
    def __div__(self, other):
        q = other.p1^2-2*other.p2^2
        return qfp2((self.p1*other.p1-2*self.p2*other.p2)/q, (self.p2*other.p1-self.p1*other.p2)/q)
    def sqrt(self):
        if is_square(self.p1):
            return qfp2(sqrt(self.p1))
        x = sqrt((self.p1 - sqrt(self.p1*self.p1 - 2*self.p2*self.p2))/2)
        y = sqrt(self.p1/2) if x == 0 else self.p2/(2*x)
        return qfp2(x, y)
    def norm(self):
        return self.p1*self.p1-2*self.p2*self.p2

def is_qfp2square(self):
    a,b = self.p1, self.p2
    def is_int(x):
        return (x - floor(x) == 0)
    if b == 0:
        return is_square(a)
    d = a*a-2*b*b
    if not is_square(d):
        return False
    c = sqrt(d)
    x = sqrt((a-c)/2)
    if not is_square(x):
        x = sqrt((a+c)/2)
        if not is_square(x):
            return False
        y = b/x
    else:
        y = 1
    if is_int(y) and ((is_int(c) and is_int(x) and x!= 0) or (x == 0)):
        return True
    else:
        return False

def GCD_qfp2(self, other):
    j,k = self, other
    if j.norm() < k.norm():
        j,k = k,j
    #print "j", j, "k:", k
    while true:
        q = j/k
        qq = q.ni()
        #print "j/k", q
        k1 = k*qq
        remainder = j - k1
        #print "qq", qq, "remainder", remainder, "k:", k, "k1:", k1
        if remainder == qfp2(0,0):
            break
        j,k = k,remainder
        #print "j", j, "k:", k
    return qfp2(abs(k.p1), abs(k.p2))


fl = []
sc = []
success_count = failure_count = 0
zero = qfp2(0)
one = qfp2(1)

for m, n in enumerate(failuresPlus2):
    #print "{0}    {1}".format(m[1], uvMinus2[m[0]][0])
    a = qfp2(n, 0)
    b = qfp2(uvPlus2[m][0], 1)
    c = GCD_qfp2(a,b)
    print "\nN:", n
    #print "a,b,c",a,b,c
    if c.norm() < 0:
        c = c * qfp2(1,1)
        print "a,b,c",a,b,c

    x1 = float(sqrt((c.p1 + float(sqrt(n)))/2))
    x2 = float(abs(c.p2/(2*x1)))
    #print "x1:", x1, "x2:", x2
    sq = qfp2(x1,x2)
    rounded = sq.ni()
    #print "{0}^2 - {1}*{2}^2 = -2".format(b.p1, m[1], uvMinus2[m[0]][1])
    #print "HCF(({0},0),({1},{2})) = ({3},{4})      sq: {5}  rounded: {6} rounded^2: {7}\n".format(a, b.p1, b.p2, c.p1, c.p2, sq, rounded, rounded*rounded)
    A = [0, rounded]
    B = [0, rounded]
    C = [one, c - rounded*rounded]
    P = [one, rounded]
    Q = [zero, one]
    i = 1
    #print "A: {0:<8}   B: {1:<8}   C: {2:<8}   P: {3:<10}  i: {4}".format(A[-1], B[-1], C[-1], P[-1], i)
    s = set()
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
        a1, a2, a3, a4 = A[-1].qftuple(), B[-1].qftuple(), C[-1].qftuple(), P[-1].qftuple()
        a5 = (a1[0], a1[1], a2[0], a2[1], a3[0], a3[1], a4[0], a4[1], )
        if a5 in s:
            str = "i: {0:<3} Repeated A, B, C, P".format(i)
            fl.append((n, str))
            failure_count += 1
            break
        else:
            s.add(a5)
        ct = -C[-1] if i%2 > 0 else C[-1] 
        if is_qfp2square(ct):
            cs = sqrt(ct)
            str = "i: {0:<3} sqrt{1} = {2}".format(i, ct, cs)
            print str
            factor = []
            t1 = P[-1] + cs
            if t1 != zero:
                f1 = GCD_qfp2(t1, c)
                factor.append(f1.norm())
            t2 = P[-1] - cs
            if t2 != zero:
                f2 = GCD_qfp2(t2, c)
                factor.append(f2.norm())
#            if len(factor) == 2:
#                print "{0} - factors {1}".format(a.p1, factor)
            if len(factor) == 2 and abs(factor[0])*abs(factor[1]) == a.p1:
                print "{0} = {1}*{2}".format(a.p1, abs(factor[0]), abs(factor[1]))
                sc.append(n)
                success_count += 1
                break
            if C[-1] == qfp2(1,0):
                fl.append((n, str))
                failure_count += 1
                print "\n"
                break
print "\nsuccess_count: {0}    failure_count: {1}   failure percentage: {2:.2f}%\n".format(success_count, failure_count, float((failure_count * 100)/(success_count+failure_count)))

now_cpu = time.clock(); now_time = time.time()
print "Time to end = {0:.2f}".format(now_time-mid_time)

#print "Success list:"
#for s in sc:
#    print s
#print "\nFailure list:"
#for s in fl:
#    print s
︡af5c8d07-573f-4b42-9b18-90eaae879be7︡{"stdout":"SageMath version 7.3, Release Date: 2016-08-04\n\n"}︡{"stdout":"N:    1001   count:     91   failuresEvenLen:     9(9.89%)   failuresOddLen:     2(2.20%)   Success: 87.91%   time: 0.18"}︡{"stdout":"\nN:    2009   count:    221   failuresEvenLen:    18(8.14%)   failuresOddLen:     4(1.81%)   Success: 90.05%   time: 0.41"}︡{"stdout":"\nN:    5017   count:    649   failuresEvenLen:    47(7.24%)   failuresOddLen:    15(2.31%)   Success: 90.45%   time: 0.99"}︡{"stdout":"\nN:   10001   count:   1413   failuresEvenLen:    90(6.37%)   failuresOddLen:    28(1.98%)   Success: 91.65%   time: 1.73"}︡{"stdout":"\nN:   20003   count:   3036   failuresEvenLen:   173(5.70%)   failuresOddLen:    70(2.31%)   Success: 92.00%   time: 4.05"}︡{"stdout":"\nEven Failures(173): \n  [119, 287, 343, 527, 623, 731, 779, 803, 959, 1003, 1127, 1331, 1343, 1411, 1519, 1679, 1687, 1691, 2023, 2123, 2147, 2263, 2599, 2603, 2651, 2747, 2807, 2863, 3239, 3247, 3403, 3479, 3587, 3647, 3667, 3707, 3791, 3827, 3971, 4151, 4223, 4247, 4319, 4487, 4499, 4579, 4763, 5027, 5327, 5339, 5371, 5383, 5543, 5627, 6167, 6191, 6239, 6319, 6347, 6443, 6499, 6503, 6559, 6611, 6727, 6859, 6887, 7223, 7339, 7343, 7403, 7567, 7571, 7619, 8023, 8227, 8279, 8459, 8483, 8519, 8531, 8651, 8687, 8899, 8927, 9023, 9079, 9247, 9407, 9799, 10199, 10291, 10307, 10327, 10423, 10447, 10727, 10747, 10823, 10951, 10999, 11023, 11099, 11207, 11263, 11303, 11327, 11543, 11711, 11723, 11879, 12167, 12179, 12319, 12431, 12587, 12767, 13067, 13123, 13223, 13379, 13391, 13439, 13739, 13823, 13943, 14003, 14191, 14267, 14507, 14603, 14611, 14903, 14959, 15127, 15479, 15623, 15967, 16019, 16199, 16379, 16439, 16571, 16639, 16643, 16847, 17063, 17083, 17111, 17311, 17363, 17639, 17687, 17867, 18071, 18227, 18599, 18631, 18767, 18779, 18983, 18991, 19103, 19307, 19439, 19547, 19567, 19723, 19771, 19871, 19879, 19883, 19967]   \nOdd Failures(70): \n  [533, 697, 1313, 1853, 2197, 2249, 2813, 3029, 3281, 3293, 3341, 3653, 3973, 4469, 4913, 5213, 5353, 5629, 5713, 5837, 5933, 6253, 6893, 8321, 8653, 9197, 9773, 9953, 10397, 11029, 11141, 11401, 11453, 11509, 11629, 11773, 11849, 12053, 12389, 12461, 12557, 12629, 12773, 13549, 13801, 14089, 14453, 15133, 15509, 15529, 15613, 15977, 16109, 16133, 16601, 16609, 16733, 17177, 17261, 17693, 17797, 17849, 17893, 18317, 18409, 18989, 19109, 19409, 19633, 19637]\n"}︡{"stdout":"failuresPlus2(103):\n [119, 287, 343, 527, 623, 959, 1127, 1343, 1519, 1679, 1687, 2023, 2263, 2599, 2807, 2863, 3239, 3247, 3479, 3647, 3791, 4151, 4223, 4247, 4319, 4487, 5327, 5383, 5543, 6167, 6191, 6239, 6319, 6503, 6559, 6727, 6887, 7223, 7343, 7567, 8023, 8279, 8519, 8687, 8927, 9023, 9079, 9247, 9407, 9799, 10199, 10327, 10423, 10447, 10727, 10823, 10951, 10999, 11023, 11207, 11263, 11303, 11327, 11543, 11711, 11879, 12167, 12319, 12431, 12767, 13223, 13391, 13439, 13823, 13943, 14191, 14903, 14959, 15127, 15479, 15623, 15967, 16199, 16439, 16639, 16847, 17063, 17111, 17311, 17639, 17687, 18071, 18599, 18631, 18767, 18983, 18991, 19103, 19439, 19567, 19871, 19879, 19967]\n"}︡{"stdout":"uvPlus2:\n[(11, 1), (17, 1), (11427, 617), (23, 1), (25, 1), (31, 1), (235, 7), (623, 17), (39, 1), (41, 1), (5886885, 143327), (45, 1), (333, 7), (51, 1), (53, 1), (19333527, 361327), (339709, 5969), (57, 1), (59, 1), (28927, 479), (431, 7), (451, 7), (65, 1), (4627, 71), (5849, 89), (67, 1), (73, 1), (7557, 103), (8413, 113), (22067, 281), (131281409, 1668487), (79, 1), (61408551, 772511), (4538405, 56279), (81, 1), (737201475, 8988257), (83, 1), (85, 1), (533257, 6223), (87, 1), (627, 7), (91, 1), (17629, 191), (4567, 49), (58063745, 614543), (95, 1), (14738061531, 154675439), (376628642558367, 3916633337039), (97, 1), (99, 1), (101, 1), (12093, 119), (365187, 3577), (1108964158695, 10849799503), (725, 7), (8177773, 78607), (1779, 17), (1325739, 12641), (105, 1), (20914110437, 197558041), (3124869580137, 29444545447), (75803, 713), (745, 7), (7843, 73), (2489, 23), (109, 1), (15191131279719638885, 137720427724123513), (111, 1), (24554054338961, 220226805607), (113, 1), (115, 1), (1483871, 12823), (372838925159, 3216160919), (823, 7), (127645, 1081), (15129, 127), (20387, 167), (40239, 329), (123, 1), (5101, 41), (125, 1), (8351025611469948705, 66088843624856687), (38213261, 300241), (4627139, 36089), (129, 1), (66637393, 513401), (11661323, 89273), (389121675379, 2974731193), (921, 7), (145117280701, 1092652721), (133, 1), (941, 7), (16229, 119), (19728190029, 144533713), (137, 1), (11216435, 81409), (182871, 1327), (119452903, 864263), (138169, 991), (7254705, 51863), (13044187034791, 92535239993), (141, 1), (21337, 151)]\n"}︡{"stdout":"Time to split = 6.32\n"}︡{"stdout":"\nN: 119\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 287\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 343\ni: 4   sqrt(1, 0) = (1, 0)\n\n\n\nN: 527\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 623\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 959\ni: 6   sqrt(1, 0) = (1, 0)\n\n\n\nN: 1127\na,b,c (1127, 0) (235, 1) (53, 29)\n\nN: 1343\ni: 6   sqrt(3, 2) = (1, 1)\n1343 = 79*17\n\nN: 1519\ni: 8   sqrt(1, 0) = (1, 0)\n\n\n\nN: 1679\ni: 16  sqrt(3, -2) = (1, -1)\n1679 = 23*73\n\nN: 1687\ni: 5   sqrt(1, 0) = (1, 0)\n1687 = 7*241\n\nN: 2023\ni: 2   sqrt(1, 0) = (1, 0)\n2023 = 289*7\n\nN: 2263\na,b,c (2263, 0) (333, 1) (75, 41)\ni: 23  sqrt(3, 2) = (1, 1)"}︡{"stdout":"\n2263 = 31*73\n\nN: 2599\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 2807\ni: 2   sqrt(1, 0) = (1, 0)\n2807 = 401*7\n\nN: 2863\ni: 6   sqrt(1, 0) = (1, 0)\n\n\n\nN: 3239\ni: 2   sqrt(3, -2) = (1, -1)\n3239 = 79*41\n\nN: 3247\ni: 9   sqrt(1, 0) = (1, 0)\ni: 18  sqrt(1, 0) = (1, 0)\n\n\n\nN: 3479\ni: 4   sqrt(1, 0) = (1, 0)\n\n\n\nN: 3647\na,b,c (3647, 0) (28927, 1) (103, 59)\ni: 6   sqrt(1, 0) = (1, 0)\n3647 = 521*7\n\nN: 3791\na,b,c (3791, 0) (431, 1) (97, 53)\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 4151\na,b,c (4151, 0) (451, 1) (101, 55)\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 4223\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 4247\ni: 14  sqrt(1, 0) = (1, 0)\n\n\n\nN: 4319\na,b,c (4319, 0) (5849, 1) (113, 65)\ni: 6   sqrt(3, 2) = (1, 1)\n4319 = 7*617\n\nN: 4487\ni: 3   sqrt(9, -4) = (3, 0)\ni: 5   sqrt(9, -4) = (3, 0)\ni: 8   sqrt(1, 0) = (1, 0)\n4487 = 7*641\n\nN: 5327\ni: 10  sqrt(1, 0) = (1, 0)\n\n\n\nN: 5383\na,b,c (5383, 0) (7557, 1) (109, 57)\ni: 8   sqrt(1, 0) = (1, 0)\n\n\n\nN: 5543\ni: 11  sqrt(3, -2) = (1, -1)\n\nN: 6167\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 6191\n\nN: 6239\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 6319\ni: 22  sqrt(3, 2) = (1, 1)\ni: 38  sqrt(3, 2) = (1, 1)\n\nN: 6503\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 6559\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 6727\ni: 4   sqrt(9, 0) = (3, 0)\ni: 10  sqrt(1, 0) = (1, 0)\n\n\n\nN: 6887\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 7223\ni: 6   sqrt(1, 0) = (1, 0)\n\n\n\nN: 7343\ni: 4   sqrt(1, 0) = (1, 0)\n7343 = 1049*7"}︡{"stdout":"\n\nN: 7567\n\nN: 8023\na,b,c (8023, 0) (627, 1) (141, 77)\ni: 5   sqrt(3, 2) = (1, 1)\n\nN: 8279\n\nN: 8519\ni: 6   sqrt(1, 0) = (1, 0)\n\n\n\nN: 8687\ni: 6   sqrt(1, 0) = (1, 0)\n8687 = 7*1241\n\nN: 8927\ni: 12  sqrt(3, 2) = (1, 1)\ni: 18  sqrt(9, 4) = (3, 0)\ni: 20  sqrt(9, 4) = (3, 0)\ni: 28  sqrt(3, -2) = (1, -1)\ni: 42  sqrt(3, 2) = (1, 1)\n\nN: 9023\n\nN: 9079\na,b,c (9079, 0) (14738061531, 1) (149, 81)\n\nN:"}︡{"stdout":" 9247\ni: 8   sqrt(1, 0) = (1, 0)\n\n\n\nN: 9407\n\nN: 9799\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 10199\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 10327\ni: 4   sqrt(9, -4) = (3, 0)\ni: 10  sqrt(3, -2) = (1, -1)\ni: 18  sqrt(3, -2) = (1, -1)\n10327 = 23*449\n\nN: 10423\ni: 17  sqrt(3, -2) = (1, -1)"}︡{"stdout":"\n10423 = 7*1489\n\nN: 10447\ni: 10  sqrt(1, 0) = (1, 0)\n\n\n\nN: 10727\na,b,c (10727, 0) (725, 1) (163, 89)\ni: 10  sqrt(1, 0) = (1, 0)\n\n\n\nN: 10823\na,b,c (10823, 0) (8177773, 1) (179, 103)\ni: 12  sqrt(3, -2) = (1, -1)\n10823 = 137*79\n\nN: 10951\ni: 7   sqrt(3, 2) = (1, 1)\ni: 23  sqrt(9, 4) = (3, 0)\ni: 30  sqrt(3, -2) = (1, -1)\n10951 = 47*233\n\nN: 10999\ni: 9   sqrt(3, -2) = (1, -1)\ni: 15  sqrt(9, 0) = (3, 0)\ni: 18  sqrt(3, 2) = (1, 1)\n10999 = 647*17\n\nN: 11023\ni: 14  sqrt(3, -2) = (1, -1)\n11023 = 73*151\n\nN: 11207\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 11263\na,b,c (11263, 0) (3124869580137, 1) (169, 93)\ni: 8   sqrt(1, 0) = (1, 0)\n11263 = 1609*7\n\nN: 11303\ni: 12  sqrt(1, 0) = (1, 0)\n11303 = 127*89\n\nN: 11327\na,b,c (11327, 0) (745, 1) (167, 91)\n\nN: 11543\na,b,c (11543, 0) (7843, 1) (181, 103)\n\nN: 11711\ni: 3   sqrt(9, 4) = (3, 0)\ni: 9   sqrt(9, 4) = (3, 0)\ni: 12  sqrt(1, 0) = (1, 0)\n11711 = 239*49\n\nN: 11879\ni: 4   sqrt(1, 0) = (1, 0)\n11879 = 7*1697\n\nN: 12167\na,b,c (12167, 0) (15191131279719638885, 1) (157, 79)\ni: 13  sqrt(9, 4) = (3, 0)\ni: 22  sqrt(3, -2) = (1, -1)\ni: 36  sqrt(9, 4) = (3, 0)\ni: 45  sqrt(3, -2) = (1, -1)\ni: 59  sqrt(9, 4) = (3, 0)\ni: 68  sqrt(3, -2) = (1, -1)\ni: 82  sqrt(9, 4) = (3, 0)"}︡{"stdout":"\ni: 91  sqrt(3, -2) = (1, -1)\n\nN: 12319\ni: 11  sqrt(9, 4) = (3, 0)\ni: 13  sqrt(9, 0) = (3, 0)\ni: 23  sqrt(3, -2) = (1, -1)\ni: 47  sqrt(3, -2) = (1, -1)\n\nN: 12431\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 12767\n\nN: 13223\ni: 9   sqrt(3, -2) = (1, -1)\n\nN: 13391\ni: 6   sqrt(1, 0) = (1, 0)\n\n\n\nN: 13439\na,b,c (13439, 0) (372838925159, 1) (167, 85)\ni: 2   sqrt(3, -2) = (1, -1)\n13439 = 151*89\n\nN: 13823\na,b,c (13823, 0) (823, 1) (185, 101)\n\nN: 13943\na,b,c (13943, 0) (127645, 1) (181, 97)\ni: 16  sqrt(9, -4) = (3, 0)\ni: 21  sqrt(3, -2) = (1, -1)\n13943 = 191*73\n\nN: 14191\n\nN: 14903\ni: 2   sqrt(1, 0) = (1, 0)\n14903 = 2129*7\n\nN: 14959\ni: 4   sqrt(1, 0) = (1, 0)\n14959 = 7*2137\n\nN: 15127\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 15479\na,b,c (15479, 0) (5101, 1) (179, 91)\ni: 24  sqrt(3, 2) = (1, 1)"}︡{"stdout":"\n15479 = 673*23\n\nN: 15623\n\nN: 15967\ni: 5   sqrt(9, 0) = (3, 0)\n\nN: 16199\na,b,c (16199, 0) (38213261, 1) (181, 91)\ni: 4   sqrt(3, -2) = (1, -1)\n16199 = 97*167\n\nN: 16439\ni: 6   sqrt(1, 0) = (1, 0)\n\n\n\nN: 16639\ni: 8   sqrt(1, 0) = (1, 0)\n\n\n\nN: 16847\na,b,c (16847, 0) (66637393, 1) (193, 101)\ni: 5   sqrt(9, -4) = (3, 0)\n\nN: 17063\ni: 2   sqrt(1, 0) = (1, 0)\n\n\n\nN: 17111\ni: 2   sqrt(9, -4) = (3, 0)\n\nN: 17311\na,b,c (17311, 0) (921, 1) (207, 113)\ni: 6   sqrt(1, 0) = (1, 0)\n\n\n\nN: 17639\ni: 22  sqrt(1, 0) = (1, 0)\n17639 = 31*569\n\nN: 17687\ni: 4   sqrt(1, 0) = (1, 0)\n17687 = 769*23\n\nN: 18071\na,b,c (18071, 0) (941, 1) (211, 115)\ni: 9   sqrt(3, -2) = (1, -1)\n18071 = 1063*17\n\nN: 18599\ni: 4   sqrt(1, 0) = (1, 0)\n18599 = 7*2657\n\nN: 18631\n\nN:"}︡{"stdout":" 18767\ni: 36  sqrt(1, 0) = (1, 0)\n\n\n\nN: 18983\na,b,c (18983, 0) (11216435, 1) (211, 113)\ni: 8   sqrt(1, 0) = (1, 0)\n\n\n\nN: 18991\n\nN:"}︡{"stdout":" 19103\ni: 10  sqrt(1, 0) = (1, 0)\n19103 = 7*2729\n\nN: 19439\na,b,c (19439, 0) (138169, 1) (241, 139)\ni: 4   sqrt(1, 0) = (1, 0)\n19439 = 7*2777\n\nN: 19567\ni: 7   sqrt(9, 0) = (3, 0)\n19567 = 17*1151\n\nN: 19871\ni: 6   sqrt(9, 4) = (3, 0)\ni: 12  sqrt(3, 2) = (1, 1)\ni: 18  sqrt(9, 4) = (3, 0)\ni: 25  sqrt(1, 0) = (1, 0)\ni: 31  sqrt(9, 4) = (3, 0)\ni: 37  sqrt(3, 2) = (1, 1)\ni: 43  sqrt(9, 4) = (3, 0)\ni: 50  sqrt(1, 0) = (1, 0)\n\n\n\nN: 19879\ni: 22  sqrt(1, 0) = (1, 0)"}︡{"stdout":"\n19879 = 193*103\n\nN: 19967\na,b,c (19967, 0) (21337, 1) (233, 131)\ni: 38  sqrt(3, -2) = (1, -1)\n19967 = 487*41\n"}︡{"stdout":"\nsuccess_count: 37    failure_count: 66   failure percentage: 64.08%\n\n"}︡{"stdout":"Time to end = 6.93\n"}︡{"done":true}︡









