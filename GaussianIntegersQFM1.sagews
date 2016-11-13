︠a055a12b-5dff-4d03-98da-5ca910736394s︠
if '__SAGEWS__' in globals():
    sage_server.MAX_OUTPUT = 100000000
    sage_server.MAX_OUTPUT_MESSAGES = 100000
import time
from abc import ABCMeta, abstractmethod

foundCount = 0
failureList =[]

FailedMod4plus1List = [ (2, 23),  (16, 21),  (32, 17),  (2, 43),  (32, 35),  (2, 53),  (2, 55),  (16, 55),  (22, 53),  (46, 35),  (58, 17),  (2, 63),  (62, 25),  (58, 43),  (72, 13),  (2, 75),  (72, 23),  (46, 61),  (2, 77),  (18, 77), (2, 83), (20,89),  (2, 93),  (94, 19),  (98, 13),  (88, 47),  (14, 101),  (2, 105),  (70, 79),  (76, 75),  (2, 107),  (30, 103),  (102, 35),  (18, 107),  (100, 43),  (38, 103),  (110, 17),  (110, 19),  (26, 109),  (110, 23),  (2, 113),  (70, 93),  (24, 115),  (92, 75),  (62, 103),  (2, 123),  (70, 103),  (48, 115),  (122, 27),  (76, 101),  (22, 125),  (2, 127),  (124, 35),  (128, 15),  (118, 53),  (4, 131),  (10, 131),  (2, 133),  (106, 81),  (80, 107),  (118, 63),  (34, 131),  (128, 45),  (110, 83),  (122, 65),  (128, 55),  (128, 57),  (74, 119)]

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

class qfm1(qf):
    """Quadratic Field Class (-1)"""
    def ni(self):
        return qfm1(floor(self.p1+0.5), floor(self.p2+0.5))
    def __neg__(self):
        return qfm1(-self.p1, -self.p2)
    def __sub__(self, other):
        return qfm1(self.p1-other.p1, self.p2-other.p2)
    def __add__(self, other):
        return qfm1(self.p1+other.p1, self.p2+other.p2)
    def __mul__(self, other):
        return qfm1(self.p1*other.p1-self.p2*other.p2, self.p2*other.p1+self.p1*other.p2)
    def __div__(self, other):
        q = other.p1^2+other.p2^2
        return qfm1((self.p1*other.p1+self.p2*other.p2)/q, (self.p2*other.p1-self.p1*other.p2)/q)
    def sqrt(self):
        x = sqrt((self.p1 + sqrt(self.p1*self.p1 + self.p2*self.p2))/2)
        y = sqrt((-self.p1 + sqrt(self.p1*self.p1 + self.p2*self.p2))/2)
        if self.p2 < 0:
            y *= -1
        return qfm1(x, y)
    def norm(self):
        return self.p1*self.p1+self.p2*self.p2
    def ni(z):
        return qfm1(floor(z.p1+0.5), floor(z.p2+0.5))


def GCDComplex(j,k):
    j1, k1 = j, k
    #print "j1", j1, "k1", k1
    while True:
        q = j/k
        #print ("q", q,"j",j,"k",k)
        qq = qfm1(round(q.p1),round(q.p2))
        remainder = j-k*qq
        if remainder == qfm1(0, 0):
            k = qfm1(abs(k.p1),abs(k.p2))
            #print 'GCD({0},{1}) = {2}'.format(j1, k1,k)
            break
        j,k = k,remainder
    return k

def fp_complex(n):
    #print n
    global foundCount, failureList

    #def pr(z):
    #    return str((z.p1, z.p2))

    def a(x):
        """ Determines whether both real and imaginary parts of complex number x
        are integers - or very close to integers"""
        rf = abs(abs(x.p1) - floor(abs(x.p1)+0.5))
        cf = abs(abs(x.p2) - floor(abs(x.p2)+0.5))
        if (rf < 0.0000001) and (cf < 0.0000001):
            return (True, qfm1(floor(x.p1+0.5), floor(x.p2+0.5)))
        else:
            return (False, x)
    tmp = sqrt(n)
    sq = tmp.ni()
    #print sq
    A = [0, sq]
    B = [0, sq]
    C = [qfm1(1,0), n - sq*sq]
    P = [qfm1(1,0),sq]
    Q = [qfm1(0,0), qfm1(1,0)]
    i = 1
    saved_i = -10
    while saved_i != i-2:
        if C[i] == qfm1(1,0) and saved_i == -10:
            saved_i = i
        k = i-1
        j = i
        i += 1
        A.append(((tmp+B[j])/C[j]).ni())
        B.append(A[i]*C[j]-B[j])
        C.append(C[k]+A[i]*(B[j]-B[i]))
        P.append(P[k]+A[i]*P[j])
        Q.append(Q[k]+A[i]*Q[j])
        a1 = sqrt(C[i]).p1
        a2 = N(a1)
        b1 = sqrt(C[i]).p2
        b2 = N(b1)
        f = a(qfm1(a2,b2))
        #print"i: {0:2}   A: {1:12} B: {2:12} C: {3:12} P: {4:30}  Q: {5:30}  {6}".format(i, pr(A[i]), pr(B[i]), pr(C[i]), pr(P[i]), pr(Q[i]), pr(f[1]))
        if f[0]:
            #print"i: {0:2}   A: {1:12} B: {2:12} C: {3:12} P: {4:30}  Q: {5:30}  {6}".format(i, pr(A[i]), pr(B[i]), pr(C[i]), pr(P[i]), pr(Q[i]), pr(f[1]))
            ff = f[1]
            if i % 2 == 1:
                ff = qfm1(-f[1].p2, f[1].p1)
            nf = qfm1(n.norm(),0)
            #print "P[-1]", P[-1], "ff", ff, "nf", nf
            f1 = GCDComplex(P[-1] + ff, nf)
            f2 = GCDComplex(P[-1] - ff, nf)
            foundFactor = false
            for f in [f1, f2]:
                #print "f.norm()", f.norm(), "nf", nf
                if f.norm() != 1 and f.norm() < n.norm():
                    print "Found factor: {0} of {1}".format(norm(f), nf)
                    foundFactor = true
            if foundFactor:
                foundCount += 1
                return True
    if saved_i == i-2:
        print 'failure:', n.norm()
        oddEven = ('even', 'odd')[int(n.p1)%2 == 0]
        failureList.append((n.norm(), oddEven))
        return False

start_cpu = time.clock(); start_time = time.time()

for f in FailedMod4plus1List:
    n = f[0]^2+f[1]^2
    print 'n:', n
    if not fp_complex(qfm1(f[1],f[0])):
        if fp_complex(qfm1(f[0],f[1])):
            failureList = failureList[:-1]
    print


#print "Found a factor in {0:.2f}% of the cases".format(float(100.0*foundCount/len(FailedMod4plus1List)))
print "Failure List:", failureList

now_cpu = time.clock(); now_time = time.time()
print "Elapsed time {0:.2f}".format(now_time-start_time)
︡ca1fc98d-9648-4ffc-a212-a4e5b6db3db6︡{"stdout":"n: 533\nFound factor: 13 of (533, 0)\nFound factor: 41 of (533, 0)\n\nn: 697\nFound factor: 41 of (697, 0)\nFound factor: 17 of (697, 0)\n\nn: 1313\nFound factor: 101 of (1313, 0)\n\nn: 1853\nFound factor: 109 of (1853, 0)\nFound factor: 17 of (1853, 0)\n\nn: 2249\nFound factor: 173 of (2249, 0)\nFound factor: 13 of (2249, 0)\n\nn: 2813\nFound factor: 97 of (2813, 0)\nFound factor: 29 of (2813, 0)\n\nn: 3029\nFound factor: 13 of (3029, 0)\nFound factor: 233 of (3029, 0)\n\nn: 3281\nfailure: 3281\nfailure:"}︡{"stdout":" 3281\n\nn: 3293\nFound factor: 89 of (3293, 0)\nFound factor: 37 of (3293, 0)\n\nn: 3341\nFound factor: 13 of (3341, 0)\nFound factor: 257 of (3341, 0)\n\nn: 3653\nFound factor: 13 of (3653, 0)\nFound factor: 281 of (3653, 0)\n\nn: 3973\nFound factor: 29 of (3973, 0)\nFound factor: 137 of (3973, 0)\n\nn: 4469\nFound factor: 41 of (4469, 0)\nFound factor: 109 of (4469, 0)\n\nn: 5213\nFound factor: 13 of (5213, 0)\nFound factor: 401 of (5213, 0)\n\nn: 5353\nFound factor: 53 of (5353, 0)\nFound factor: 101 of (5353, 0)\n\nn: 5629\nFound factor: 13 of (5629, 0)\nFound factor: 433 of (5629, 0)\n\nn: 5713\nFound factor: 29 of (5713, 0)\nFound factor: 197 of (5713, 0)\n\nn: 5837\nFound factor: 449 of (5837, 0)\nFound factor: 13 of (5837, 0)\n\nn: 5933\nFound factor: 349 of (5933, 0)"}︡{"stdout":"\nFound factor: 17 of (5933, 0)\n\nn: 6253\nFound factor: 169 of (6253, 0)\nFound factor: 37 of (6253, 0)\n\nn: 6893\nfailure: 6893\nFound factor: 113 of (6893, 0)\nFound factor: 61 of (6893, 0)\n\nn: 8321\nFound factor: 157 of (8321, 0)\nFound factor: 53 of (8321, 0)\n\nn: 8653\nFound factor: 509 of (8653, 0)\nFound factor: 17 of (8653, 0)\n\nn: 9197\nFound factor: 17 of (9197, 0)\nFound factor: 541 of (9197, 0)\n\nn: 9773\nfailure:"}︡{"stdout":" 9773\nFound factor: 29 of (9773, 0)\nFound factor: 337 of (9773, 0)\n\nn: 9953\nFound factor: 37 of (9953, 0)\n\nn: 10397\nFound factor: 281 of (10397, 0)\nFound factor: 37 of (10397, 0)\n\nn: 11029\nfailure:"}︡{"stdout":" 11029\nFound factor: 41 of (11029, 0)\nFound factor: 269 of (11029, 0)\n\nn: 11141\nFound factor: 857 of (11141, 0)\nFound factor: 13 of (11141, 0)\n\nn: 11401\nFound factor: 877 of (11401, 0)\nFound factor: 13 of (11401, 0)\n\nn: 11453\nFound factor: 13 of (11453, 0)\nFound factor: 881 of (11453, 0)\n\nn: 11509\nFound factor: 17 of (11509, 0)\nFound factor: 677 of (11509, 0)\n\nn: 11629\nFound factor: 29 of (11629, 0)\nFound factor: 401 of (11629, 0)\n\nn: 11773\nFound factor: 193 of (11773, 0)\nFound factor: 61 of (11773, 0)\n\nn: 11849\nFound factor: 41 of (11849, 0)"}︡{"stdout":"\nFound factor: 289 of (11849, 0)\n\nn: 12053\nFound factor: 17 of (12053, 0)\nFound factor: 709 of (12053, 0)\n\nn: 12389\nfailure: 12389\nFound factor: 13 of (12389, 0)"}︡{"stdout":"\n\nn: 12461\nFound factor: 733 of (12461, 0)\nFound factor: 17 of (12461, 0)\n\nn: 12557\nFound factor: 433 of (12557, 0)\nFound factor: 29 of (12557, 0)\n\nn: 12629\nFound factor: 73 of (12629, 0)\nFound factor: 173 of (12629, 0)\n\nn: 12773\nFound factor: 241 of (12773, 0)\nFound factor: 53 of (12773, 0)\n\nn: 13549\nFound factor: 289 of (13549, 0)"}︡{"stdout":"\nFound factor: 797 of (13549, 0)\n\nn: 13801\nFound factor: 37 of (13801, 0)\nFound factor: 373 of (13801, 0)\n\nn: 14089\nFound factor: 193 of (14089, 0)\nFound factor: 73 of (14089, 0)\n\nn: 14453\nFound factor: 97 of (14453, 0)"}︡{"stdout":"\nFound factor: 149 of (14453, 0)\n\nn: 15133\nfailure: 15133\nFound factor: 37 of (15133, 0)\nFound factor: 409 of (15133, 0)\n\nn: 15509\nFound factor: 1193 of (15509, 0)\nFound factor: 13 of (15509, 0)\n\nn: 15529\nFound factor: 53 of (15529, 0)\nFound factor: 293 of (15529, 0)\n\nn: 15613\nFound factor: 1201 of (15613, 0)"}︡{"stdout":"\nFound factor: 13 of (15613, 0)\n\nn: 15977\nFound factor: 13 of (15977, 0)\n\nn: 16109\nFound factor: 181 of (16109, 0)\nFound factor: 89 of (16109, 0)\n\nn: 16133\nFound factor: 949 of (16133, 0)\nFound factor: 221 of (16133, 0)\n\nn: 16601\nFound factor: 1277 of (16601, 0)\nFound factor: 13 of (16601, 0)\n\nn: 16609\nFound factor: 17 of (16609, 0)"}︡{"stdout":"\nFound factor: 977 of (16609, 0)\n\nn: 16733\nfailure:"}︡{"stdout":" 16733\nFound factor: 577 of (16733, 0)\n\nn: 17177\nFound factor: 193 of (17177, 0)\nFound factor: 89 of (17177, 0)\n\nn: 17261\nFound factor: 421 of (17261, 0)\nFound factor: 41 of (17261, 0)\n\nn: 17693\nFound factor: 1361 of (17693, 0)\nFound factor: 13 of (17693, 0)\n\nn: 17797\nFound factor: 169 of (17797, 0)\nFound factor: 1369 of (17797, 0)\n\nn: 17849\nFound factor: 1373 of (17849, 0)"}︡{"stdout":"\nFound factor: 13 of (17849, 0)\n\nn: 17893\nFound factor: 617 of (17893, 0)\nFound factor: 29 of (17893, 0)\n\nn: 18317\nFound factor: 13 of (18317, 0)\nFound factor: 1409 of (18317, 0)\n\nn: 18409\nFound factor: 449 of (18409, 0)"}︡{"stdout":"\nFound factor: 41 of (18409, 0)\n\nn: 18989\nFound factor: 1117 of (18989, 0)\nFound factor: 17 of (18989, 0)\n\nn: 19109\nFound factor: 97 of (19109, 0)\nFound factor: 197 of (19109, 0)\n\nn: 19409\nFound factor: 13 of (19409, 0)\n\nn: 19633\nFound factor: 29 of (19633, 0)"}︡{"stdout":"\nFound factor: 677 of (19633, 0)\n\nn: 19637\nFound factor: 73 of (19637, 0)\nFound factor: 269 of (19637, 0)\n\n"}︡{"stdout":"Failure List: [(3281, 'even'), (3281, 'odd')]\n"}︡{"stdout":"Elapsed time 1.76\n"}︡{"done":true}︡

n(3.36932485156297, 0.593590730520460) (False, (3.36932485156297, 0.593590730520460))\n(15, 1)\n(sqrt(1/2*sqrt(226) + 15/2), sqrt(1/2*sqrt(226) - 15/2))\n(3.87513202216214, 0.129027862054884) (False, (3.87513202216214, 0.129027862054884))\n(3, 4)\n(2, 1)\n(2.00000000000000, 1.00000000000000) (True, (2, 1))\nf[1]: (2, 1) ff (-1, 2)\nnf1: 19637\nnf: (19637, 0)\nj (323270332414, 1765542754820) k (19637, 0)\nj (323270332414, 1765542754820) k (19637, 0)\nGCD((323270332414, 1765542754820),(19637, 0)) = (10, 13)\nj (323270332416, 1765542754816) k (19637, 0)\nj (323270332416, 1765542754816) k (19637, 0)\nGCD((323270332416, 1765542754816),(19637, 0)) = (8, 3)\n(-13, 5)\n(sqrt(1/2*sqrt(194) - 13/2), sqrt(1/2*sqrt(194) + 13/2))\n(0.681317942367629, 3.66935881845753) (False, (0.681317942367629, 3.66935881845753))\n(-1, -8)\n(sqrt(1/2*sqrt(65) - 1/2), sqrt(1/2*sqrt(65) + 1/2))\n(1.87912981833328, 2.12864484453120) (False, (1.87912981833328, 2.12864484453120))\n(-2, -7)\n(sqrt(1/2*sqrt(53) - 1), sqrt(1/2*sqrt(53) + 1))\n(1.62482458888345, 2.15407867652049) (False, (1.62482458888345, 2.15407867652049))\n(1, -11)\n(sqrt(1/2*sqrt(122) + 1/2), sqrt(1/2*sqrt(122) - 1/2))\n(2.45411501535556, 2.24113375517697) (False, (2.45411501535556, 2.24113375517697))\n(3, 8)\n(sqrt(1/2*sqrt(73) + 3/2), sqrt(1/2*sqrt(73) - 3/2))\n(2.40249908900269, 1.66493299344411) (False, (2.40249908900269, 1.66493299344411))\n(10, -3)\n(sqrt(1/2*sqrt(109) + 5), sqrt(1/2*sqrt(109) - 5))\n(3.19689744196702, 0.469204917339189) (False, (3.19689744196702, 0.469204917339189))\n(1, 8)\n(sqrt(1/2*sqrt(65) + 1/2), sqrt(1/2*sqrt(65) - 1/2))\n(2.12864484453120, 1.87912981833328) (False, (2.12864484453120, 1.87912981833328))\n(7, 3)\n(sqrt(1/2*sqrt(58) + 7/2), sqrt(1/2*sqrt(58) - 7/2))\n(2.70331029534753, 0.554875258893343) (False, (2.70331029534753, 0.554875258893343))\n(-2, -1)\n(sqrt(1/2*sqrt(5) - 1), sqrt(1/2*sqrt(5) + 1))\n(0.343560749722513, 1.45534669022535) (False, (0.343560749722513, 1.45534669022535))\n(7, -8)"}









