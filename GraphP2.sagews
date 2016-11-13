︠3e99d8f8-3481-49e8-8679-912a09b1278bs︠
import time
from abc import ABCMeta, abstractmethod

print version(); print

def PrimeDivisorCount(n):
    pd = prime_divisors(n)
    return len(pd)

F = (x for x in range(119,20050,2) if not is_square(x) and not is_prime(x) and (x %3 != 0) and (x%5 != 0))
#F = (x for x in range(803,20050,2) if (PrimeDivisorCount(x) > 1) and not is_square(x) and not is_prime(x) and (x %3 != 0) and (x%5 != 0))
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

time_start_part2 = time.time()
zero = qfp2(0)
def is_int(x):
    return (x - floor(x) == 0)
time_last = time.time()
for m, n in enumerate(failuresPlus2):
    f1 = zero
    f2 = zero
    x1 =ceil(sqrt(n))
    while(true):
        # x^2-2*y^2 = n
        y1 = sqrt((x1*x1-n)/2)
        #print "x1:", x1, "y1:", y1
        if is_int(y1):
            f1 = qfp2(x1,y1)
            break
        x1 += 1
    # x^2-2*y^2 = -n
    y2 = ceil(sqrt(n/2))
    while(true):
        x2 = sqrt(2*y2*y2-n)
        #print "x2:", x2, "y2:", y2
        if is_int(x2):
            f2 = qfp2(x2,y2)
            break
        y2 += 1
    if f1 != zero and f2 != zero:
        #print "n:", n, "f1:", f1, "f2:", f2
        fac1 = GCD_qfp2(f1,f2)
        fac2 = GCD_qfp2(qfp2(x1, -y1),f2)
        factor1 = abs(fac1.norm())
        factor2 = abs(fac2.norm())
#        time_now = time.time()
#        print "                                {0} = {1}*{2}     {3}".format(n, factor1, factor2, time_now - time_last)
#        time_last = time_now
        print "{0} = {1}*{2}".format(n, factor1, factor2)
print "Time for part2: {0:.2f}".format(time.time()- time_start_part2)
print "Total time: {0:.2f}".format(time.time()- start_time)
︡e376c094-1bfb-4784-9dea-d5e156d62f5a︡{"stdout":"SageMath version 7.3, Release Date: 2016-08-04\n\n"}︡{"stdout":"N:    1001   count:     91   failuresEvenLen:     9(9.89%)   failuresOddLen:     2(2.20%)   Success: 87.91%   time: 0.31"}︡{"stdout":"\nN:    2009   count:    221   failuresEvenLen:    18(8.14%)   failuresOddLen:     4(1.81%)   Success: 90.05%   time: 0.51"}︡{"stdout":"\nN:    5017   count:    649   failuresEvenLen:    47(7.24%)   failuresOddLen:    15(2.31%)   Success: 90.45%   time: 0.97"}︡{"stdout":"\nN:   10001   count:   1413   failuresEvenLen:    90(6.37%)   failuresOddLen:    28(1.98%)   Success: 91.65%   time: 2.28"}︡{"stdout":"\nN:   20003   count:   3036   failuresEvenLen:   173(5.70%)   failuresOddLen:    70(2.31%)   Success: 92.00%   time: 4.15"}︡{"stdout":"\nEven Failures(173): \n  [119, 287, 343, 527, 623, 731, 779, 803, 959, 1003, 1127, 1331, 1343, 1411, 1519, 1679, 1687, 1691, 2023, 2123, 2147, 2263, 2599, 2603, 2651, 2747, 2807, 2863, 3239, 3247, 3403, 3479, 3587, 3647, 3667, 3707, 3791, 3827, 3971, 4151, 4223, 4247, 4319, 4487, 4499, 4579, 4763, 5027, 5327, 5339, 5371, 5383, 5543, 5627, 6167, 6191, 6239, 6319, 6347, 6443, 6499, 6503, 6559, 6611, 6727, 6859, 6887, 7223, 7339, 7343, 7403, 7567, 7571, 7619, 8023, 8227, 8279, 8459, 8483, 8519, 8531, 8651, 8687, 8899, 8927, 9023, 9079, 9247, 9407, 9799, 10199, 10291, 10307, 10327, 10423, 10447, 10727, 10747, 10823, 10951, 10999, 11023, 11099, 11207, 11263, 11303, 11327, 11543, 11711, 11723, 11879, 12167, 12179, 12319, 12431, 12587, 12767, 13067, 13123, 13223, 13379, 13391, 13439, 13739, 13823, 13943, 14003, 14191, 14267, 14507, 14603, 14611, 14903, 14959, 15127, 15479, 15623, 15967, 16019, 16199, 16379, 16439, 16571, 16639, 16643, 16847, 17063, 17083, 17111, 17311, 17363, 17639, 17687, 17867, 18071, 18227, 18599, 18631, 18767, 18779, 18983, 18991, 19103, 19307, 19439, 19547, 19567, 19723, 19771, 19871, 19879, 19883, 19967]   \nOdd Failures(70): \n  [533, 697, 1313, 1853, 2197, 2249, 2813, 3029, 3281, 3293, 3341, 3653, 3973, 4469, 4913, 5213, 5353, 5629, 5713, 5837, 5933, 6253, 6893, 8321, 8653, 9197, 9773, 9953, 10397, 11029, 11141, 11401, 11453, 11509, 11629, 11773, 11849, 12053, 12389, 12461, 12557, 12629, 12773, 13549, 13801, 14089, 14453, 15133, 15509, 15529, 15613, 15977, 16109, 16133, 16601, 16609, 16733, 17177, 17261, 17693, 17797, 17849, 17893, 18317, 18409, 18989, 19109, 19409, 19633, 19637]\n"}︡{"stdout":"failuresPlus2(103):\n [119, 287, 343, 527, 623, 959, 1127, 1343, 1519, 1679, 1687, 2023, 2263, 2599, 2807, 2863, 3239, 3247, 3479, 3647, 3791, 4151, 4223, 4247, 4319, 4487, 5327, 5383, 5543, 6167, 6191, 6239, 6319, 6503, 6559, 6727, 6887, 7223, 7343, 7567, 8023, 8279, 8519, 8687, 8927, 9023, 9079, 9247, 9407, 9799, 10199, 10327, 10423, 10447, 10727, 10823, 10951, 10999, 11023, 11207, 11263, 11303, 11327, 11543, 11711, 11879, 12167, 12319, 12431, 12767, 13223, 13391, 13439, 13823, 13943, 14191, 14903, 14959, 15127, 15479, 15623, 15967, 16199, 16439, 16639, 16847, 17063, 17111, 17311, 17639, 17687, 18071, 18599, 18631, 18767, 18983, 18991, 19103, 19439, 19567, 19871, 19879, 19967]\n"}︡{"stdout":"uvPlus2:\n[(11, 1), (17, 1), (11427, 617), (23, 1), (25, 1), (31, 1), (235, 7), (623, 17), (39, 1), (41, 1), (5886885, 143327), (45, 1), (333, 7), (51, 1), (53, 1), (19333527, 361327), (339709, 5969), (57, 1), (59, 1), (28927, 479), (431, 7), (451, 7), (65, 1), (4627, 71), (5849, 89), (67, 1), (73, 1), (7557, 103), (8413, 113), (22067, 281), (131281409, 1668487), (79, 1), (61408551, 772511), (4538405, 56279), (81, 1), (737201475, 8988257), (83, 1), (85, 1), (533257, 6223), (87, 1), (627, 7), (91, 1), (17629, 191), (4567, 49), (58063745, 614543), (95, 1), (14738061531, 154675439), (376628642558367, 3916633337039), (97, 1), (99, 1), (101, 1), (12093, 119), (365187, 3577), (1108964158695, 10849799503), (725, 7), (8177773, 78607), (1779, 17), (1325739, 12641), (105, 1), (20914110437, 197558041), (3124869580137, 29444545447), (75803, 713), (745, 7), (7843, 73), (2489, 23), (109, 1), (15191131279719638885, 137720427724123513), (111, 1), (24554054338961, 220226805607), (113, 1), (115, 1), (1483871, 12823), (372838925159, 3216160919), (823, 7), (127645, 1081), (15129, 127), (20387, 167), (40239, 329), (123, 1), (5101, 41), (125, 1), (8351025611469948705, 66088843624856687), (38213261, 300241), (4627139, 36089), (129, 1), (66637393, 513401), (11661323, 89273), (389121675379, 2974731193), (921, 7), (145117280701, 1092652721), (133, 1), (941, 7), (16229, 119), (19728190029, 144533713), (137, 1), (11216435, 81409), (182871, 1327), (119452903, 864263), (138169, 991), (7254705, 51863), (13044187034791, 92535239993), (141, 1), (21337, 151)]\n"}︡{"stdout":"119 = 17*7\n287 = 7*41\n343 = 49*7\n527 = 31*17\n623 = 89*7\n959 = 7*137\n1127 = 161*7\n1343 = 79*17\n1519 = 217*7\n1679 = 23*73\n1687 = 241*7\n2023 = 7*289\n2263 = 31*73\n2599 = 113*23\n2807 = 401*7\n2863 = 409*7\n3239 = 79*41\n3247 = 191*17\n3479 = 7*497\n3647 = 521*7"}︡{"stdout":"\n3791 = 223*17\n4151 = 593*7\n4223 = 103*41\n4247 = 31*137\n4319 = 617*7\n4487 = 641*7\n5327 = 7*761\n5383 = 769*7\n5543 = 241*23"}︡{"stdout":"\n6167 = 881*7\n6191 = 151*41\n6239 = 17*367\n6319 = 71*89\n6503 = 929*7\n6559 = 937*7\n6727 = 961*7\n6887 = 97*71"}︡{"stdout":"\n7223 = 31*233\n7343 = 1049*7\n7567 = 7*1081\n8023 = 71*113\n8279 = 487*17\n8519 = 1217*7"}︡{"stdout":"\n8687 = 1241*7\n8927 = 79*113\n9023 = 1289*7\n9079 = 1297*7\n9247 = 7*1321\n9407 = 409*23\n9799 = 41*239\n10199 = 7*1457\n10327 = 449*23"}︡{"stdout":"\n10423 = 1489*7\n10447 = 337*31\n10727 = 17*631\n10823 = 79*137"}︡{"stdout":"\n10951 = 233*47\n10999 = 647*17\n11023 = 73*151\n11207 = 1601*7"}︡{"stdout":"\n11263 = 1609*7\n11303 = 127*89\n11327 = 241*47"}︡{"stdout":"\n11543 = 1649*7\n11711 = 49*239\n11879 = 1697*7\n12167 = 23*529\n12319 = 127*97"}︡{"stdout":"\n12431 = 31*401\n12767 = 17*751\n13223 = 7*1889\n13391 = 1913*7\n13439 = 89*151"}︡{"stdout":"\n13823 = 23*601\n13943 = 191*73\n14191 = 617*23\n14903 = 2129*7\n14959 = 2137*7"}︡{"stdout":"\n15127 = 2161*7\n15479 = 673*23\n15623 = 919*17\n15967 = 2281*7\n16199 = 97*167"}︡{"stdout":"\n16439 = 967*17\n16639 = 7*2377\n16847 = 17*991\n17063 = 113*151\n17111 = 241*71"}︡{"stdout":"\n17311 = 2473*7\n17639 = 569*31\n17687 = 23*769\n18071 = 17*1063\n18599 = 2657*7"}︡{"stdout":"\n18631 = 31*601\n18767 = 2681*7\n18983 = 463*41\n18991 = 2713*7\n19103 = 2729*7\n19439 = 2777*7"}︡{"stdout":"\n19567 = 1151*17\n19871 = 31*641\n19879 = 193*103\n19967 = 487*41\n"}︡{"stdout":"Time for part2: 1.84\n"}︡{"stdout":"Total time: 7.45\n"}︡{"done":true}︡









