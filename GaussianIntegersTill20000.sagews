︠321fce1a-ff9e-4a9c-b939-ad0b9ae81d8fso︠
if '__SAGEWS__' in globals():
    sage_server.MAX_OUTPUT = 100000000
    sage_server.MAX_OUTPUT_MESSAGES = 100000
import time

def format_complex(z):
    return "({0},{1})".format(int(real(z)),int(imag(z)))

def GCDComplex(j,k):
    j1, k1 = j, k
    if norm(j) < norm(k):
        j,k = k,j
    while true:
        q = j/k
        qq = round(real(q)) + I*round(imag(q))
        remainder = j-k*qq
        if remainder == 0:
            k = abs(real(k))+I*abs(imag(k))
            #print 'GCD({0},{1}) = {2}'.format(format_complex(j1), format_complex(k1), format_complex(k))
            break
        j,k = k,remainder
    return k

foundCount = 0
failureList =[]
def fp_complex(n):
    global foundCount, failureList

    def ni(z):
        return floor(real(z)+0.5)+floor(imag(z)+0.5)*I
    def pr(z):
        return str((real(z), imag(z)))


    def a(x):
        """ Determines whether both real and imaginary parts of complex number x
        are integers - or very close to integers"""
        #print 'x:',x
        rf = abs(abs(real(x)) - floor(abs(real(x))+0.5))
        cf = abs(abs(imag(x)) - floor(abs(imag(x))+0.5))
        if (rf < 0.0000001) and (cf < 0.0000001):
            return (true, floor(real(x)+0.5) + floor(imag(x)+0.5)*I)
        else:
            return (false, x)


    tmp = sqrt(n)
    sq = ni(tmp)
    A = [0, sq]
    B = [0, sq]
    C = [1, n - sq*sq]
    P = [1,sq]
    Q = [0, 1]
    i = 1
    #print"i: {0:2}   A: {1:12} B: {2:12} C: {3:12} P: {4:30}  Q: {5:30}  {6}".format(i, pr(A[i]), pr(B[i]), pr(C[i]), pr(P[i]), pr(Q[i]), pr(a(N(sqrt(C[i])))[1]))

    saved_i = -10
    while saved_i != i-2:
        if C[i] == 1 and saved_i == -10:
            saved_i = i
        k = i-1
        j = i
        i += 1
        A.append(ni((tmp+B[j])/C[j]))
        B.append(A[i]*C[j]-B[j])
        C.append(C[k]+A[i]*(B[j]-B[i]))
        P.append(P[k]+A[i]*P[j])
        Q.append(Q[k]+A[i]*Q[j])
        f = a(N(sqrt(C[i])))
        #if not f[0]:
        #    print"i: {0:2}   A: {1:12} B: {2:12} C: {3:12} P: {4:30}  Q: {5:30}  {6}".format(i, pr(A[i]), pr(B[i]), pr(C[i]), pr(P[i]), pr(Q[i]), pr(f[1]))
        if f[0]:
            #print"i: {0:2}   A: {1:12} B: {2:12} C: {3:12} P: {4:30}  Q: {5:30}  {6}".format(i, pr(A[i]), pr(B[i]), pr(C[i]), pr(P[i]), pr(Q[i]), pr(f[1]))
            #print"i: {0:2}   A: {1:12} B: {2:12} C: {3:12} P: {4:30}  Q: {5:30}  {6}".format(i, pr(A[i]), pr(B[i]), pr(C[i]), pr(P[i]), pr(Q[i]), pr(f[1]))
            #print 'P: {0}    sqrt({1}): {2}'.format(P[-1], C[-1], f[1])
            ff = f[1]
            if i % 2 == 1:
                ff = I*f[1]
            nf = norm(n)
            f1 = GCDComplex(P[-1] + ff, nf)
            f2 = GCDComplex(P[-1] - ff, nf)
            foundFactor = false
            for f in [f1, f2]:
                if norm(f) != 1 and norm(f) < nf:
                    print "Found factor: {0} of {1}".format(norm(f), nf)
                    foundFactor = true
            if foundFactor:
                foundCount += 1
                return True
    if saved_i == i-2:
        print 'failure:', n.norm()
        oddEven = ('even', 'odd')[int(real(n))%2 == 0]
        failureList.append((n.norm(), oddEven))
        return False


FailedMod4plus1List = [ (2, 23),  (16, 21),  (32, 17),  (2, 43),  (32, 35),  (2, 53),  (2, 55),  (16, 55),  (22, 53),  (46, 35),  (58, 17),  (2, 63),  (62, 25),  (58, 43),  (72, 13),  (2, 75),  (72, 23),  (46, 61),  (2, 77),  (18, 77), (2, 83), (20,89),  (2, 93),  (94, 19),  (98, 13),  (88, 47),  (14, 101),  (2, 105),  (70, 79),  (76, 75),  (2, 107),  (30, 103),  (102, 35),  (18, 107),  (100, 43),  (38, 103),  (110, 17),  (110, 19),  (26, 109),  (110, 23),  (2, 113),  (70, 93),  (24, 115),  (92, 75),  (62, 103),  (2, 123),  (70, 103),  (48, 115),  (122, 27),  (76, 101),  (22, 125),  (2, 127),  (124, 35),  (128, 15),  (118, 53),  (4, 131),  (10, 131),  (2, 133),  (106, 81),  (80, 107),  (118, 63),  (34, 131),  (128, 45),  (110, 83),  (122, 65),  (128, 55),  (128, 57),  (74, 119),  (110, 91),  (2, 143),  (104, 99),  (86, 115),  (2, 145),  (98, 107),  (120, 83),  (28, 145),  (102, 107),  (70, 131),  (76, 129),  (88, 123),  (130, 79),  (2, 153),  (128, 87),  (94, 125),  (154, 29),  (24, 155),  (2, 157),  (26, 155),  (44, 151),  (126, 95),  (106, 121),  (122, 105),  (2, 165),  (52, 157),  (58, 157),  (134, 101),  (86, 145),  (136, 105),  (138, 103),  (140, 101),  (2, 173),  (8, 173),  (118, 127),  (2, 175),  (166, 59),  (106, 141),  (110, 139),  (142, 107),  (94, 151),  (56, 169),  (178, 13),  (166, 69),  (74, 165),  (26, 181),  (182, 23),  (110, 147),  (74, 169),  (166, 81),  (2, 185),  (2, 187),  (150, 113),  (180, 59),  (52, 183),  (190, 27),  (154, 115),  (54, 185),  (188, 43),  (116, 155),  (2, 195),  (118, 157),  (166, 105),  (2, 197),  (134, 145),  (130, 149),  (166, 109),  (122, 157),  (138, 145),  (86, 181),  (200, 21),  (194, 55),  (34, 199),  (172, 107),  (78, 187),  (162, 125),  (2, 205),  (182, 95),  (190, 79),  (124, 165),  (160, 131),  (184, 101),  (30, 209),  (34, 209),  (54, 205),  (170, 127),  (206, 51),  (50, 207),  (2, 213),  (172, 127),  (214, 9),  (64, 205),  (164, 139),  (214, 29),  (174, 131),  (218, 23),  (188, 113),  (218, 25),  (182, 125),  (94, 201),  (178, 133),  (2, 223),  (12, 223),  (134, 179)]

def TupleToComplex(t):
    return t[0]+i*t[1]

start_cpu = time.clock(); start_time = time.time()

for f in FailedMod4plus1List:
    n = f[0]^2+f[1]^2
    print 'n:', n
    if n > 20000:
        break
    if not fp_complex(TupleToComplex((f[1],f[0]))):
        if fp_complex(TupleToComplex((f[0],f[1]))):
            failureList = failureList[:-1]
    print


print "Found a factor in {0:.2f}% of the cases".format(float(100.0*foundCount/len(FailedMod4plus1List)))
print "Failure List:", failureList

now_cpu = time.clock(); now_time = time.time()
print "Elapsed time {0:.2f}".format(now_time-start_time)
︡9570c6cd-435d-40cb-9a97-eb30f8f5db0a︡{"stdout":"n: 533\namended ff"}︡{"stdout":" I + 1\nFound factor: 13 of 533"}︡{"stdout":"\nFound factor: 41 of 533\n\nn: 697\namended ff"}︡{"stdout":" 2*I\nFound factor: 41 of 697\nFound factor: 17 of 697\n\nn: 1313\namended ff"}︡{"stdout":" I - 2\nFound factor: 101 of 1313"}︡{"stdout":"\n\nn: 1853\namended ff"}︡{"stdout":" I + 2\nFound factor: 109 of 1853\nFound factor: 17 of 1853\n\nn: 2249\namended ff"}︡{"stdout":" 2*I + 1\nFound factor: 173 of 2249\nFound factor: 13 of 2249\n\nn: 2813\namended ff"}︡{"stdout":" 2*I - 1\namended ff"}︡{"stdout":" I + 1\nFound factor: 97 of 2813\nFound factor: 29 of 2813\n\nn: 3029\nFound factor: 13 of 3029"}︡{"stdout":"\nFound factor: 233 of 3029\n\nn: 3281\nfailure:"}︡{"stdout":" 3281\nfailure:"}︡{"stdout":" 3281\n\nn: 3293\namended ff"}︡{"stdout":" I + 1\nFound factor: 89 of 3293\nFound factor: 37 of 3293\n\nn: 3341\nFound factor: 13 of 3341"}︡{"stdout":"\nFound factor: 257 of 3341\n\nn: 3653\nFound factor: 13 of 3653"}︡{"stdout":"\nFound factor: 281 of 3653\n\nn: 3973\namended ff"}︡{"stdout":" 2*I - 1\namended ff"}︡{"stdout":" -3\nFound factor: 29 of 3973\nFound factor: 137 of 3973\n\nn: 4469\namended ff"}︡{"stdout":" 2*I - 1\nFound factor: 41 of 4469\nFound factor: 109 of 4469\n\nn: 5213\nFound factor: 13 of 5213"}︡{"stdout":"\nFound factor: 401 of 5213\n\nn: 5353\namended ff"}︡{"stdout":" I + 2\nFound factor: 53 of 5353\nFound factor: 101 of 5353\n\nn: 5629\namended ff"}︡{"stdout":" I - 3\namended ff"}︡{"stdout":" I - 2\nFound factor: 13 of 5629\nFound factor: 433 of 5629\n\nn: 5713\nFound factor: 29 of 5713"}︡{"stdout":"\nFound factor: 197 of 5713\n\nn: 5837\namended ff"}︡{"stdout":" 3*I + 1\namended ff"}︡{"stdout":" I + 2\nFound factor: 449 of 5837\nFound factor: 13 of 5837\n\nn: 5933\namended ff"}︡{"stdout":" I - 2\namended ff"}︡{"stdout":" I - 2\nFound factor: 349 of 5933\nFound factor: 17 of 5933\n\nn: 6253\nFound factor: 169 of 6253"}︡{"stdout":"\nFound factor: 37 of 6253\n\nn: 6893\namended ff"}︡{"stdout":" I + 1\nfailure:"}︡{"stdout":" 6893\nFound factor: 113 of 6893"}︡{"stdout":"\nFound factor: 61 of 6893\n\nn: 8321\namended ff"}︡{"stdout":" -2\nFound factor: 157 of 8321\nFound factor: 53 of 8321\n\nn: 8653\namended ff"}︡{"stdout":" 2*I - 1\namended ff"}︡{"stdout":" 2*I - 1\nFound factor: 509 of 8653\nFound factor: 17 of 8653\n\nn: 9197\nFound factor: 17 of 9197"}︡{"stdout":"\nFound factor: 541 of 9197\n\nn: 9773\nfailure:"}︡{"stdout":" 9773\namended ff"}︡{"stdout":" 2*I + 1\nFound factor: 29 of 9773\nFound factor: 337 of 9773\n\nn: 9953\nFound factor: 37 of 9953"}︡{"stdout":"\n\nn: 10397\namended ff"}︡{"stdout":" I + 3\namended ff"}︡{"stdout":" I - 2\nFound factor: 281 of 10397\nFound factor: 37 of 10397\n\nn: 11029\namended ff"}︡{"stdout":" I - 1\nfailure:"}︡{"stdout":" 11029\nFound factor: 41 of 11029"}︡{"stdout":"\nFound factor: 269 of 11029\n\nn: 11141\nFound factor: 857 of 11141"}︡{"stdout":"\nFound factor: 13 of 11141\n\nn: 11401\namended ff"}︡{"stdout":" 2*I - 2\nFound factor: 877 of 11401\nFound factor: 13 of 11401\n\nn: 11453\nFound factor: 13 of 11453"}︡{"stdout":"\nFound factor: 881 of 11453\n\nn: 11509\namended ff"}︡{"stdout":" -3\nFound factor: 17 of 11509\nFound factor: 677 of 11509\n\nn: 11629\namended ff"}︡{"stdout":" -1\nFound factor: 29 of 11629\nFound factor: 401 of 11629\n\nn: 11773\namended ff"}︡{"stdout":" I - 1\nFound factor: 193 of 11773\nFound factor: 61 of 11773\n\nn: 11849\namended ff"}︡{"stdout":" 2*I - 2\nFound factor: 41 of 11849\nFound factor: 289 of 11849\n\nn: 12053\namended ff"}︡{"stdout":" I - 2\namended ff"}︡{"stdout":" I + 1\nFound factor: 17 of 12053\nFound factor: 709 of 12053\n\nn: 12389\namended ff"}︡{"stdout":" I + 1\namended ff"}︡{"stdout":" I - 1\nfailure:"}︡{"stdout":" 12389\nFound factor: 13 of 12389"}︡{"stdout":"\n\nn: 12461\namended ff"}︡{"stdout":" 2*I + 1\nFound factor: 733 of 12461"}︡{"stdout":"\nFound factor: 17 of 12461\n\nn: 12557\nFound factor: 433 of 12557"}︡{"stdout":"\nFound factor: 29 of 12557\n\nn: 12629\namended ff"}︡{"stdout":" 3*I - 2\nFound factor: 73 of 12629"}︡{"stdout":"\nFound factor: 173 of 12629\n\nn: 12773\namended ff"}︡{"stdout":" I - 3\namended ff"}︡{"stdout":" I - 3\nFound factor: 241 of 12773\nFound factor: 53 of 12773\n\nn: 13549\nFound factor: 289 of 13549"}︡{"stdout":"\nFound factor: 797 of 13549\n\nn: 13801\nFound factor: 37 of 13801"}︡{"stdout":"\nFound factor: 373 of 13801\n\nn: 14089\namended ff"}︡{"stdout":" I + 2\namended ff"}︡{"stdout":" 3*I\namended ff"}︡{"stdout":" I - 2\nFound factor: 193 of 14089\nFound factor: 73 of 14089\n\nn: 14453\namended ff"}︡{"stdout":" I - 3\nFound factor: 97 of 14453"}︡{"stdout":"\nFound factor: 149 of 14453\n\nn: 15133\namended ff"}︡{"stdout":" I + 1\nfailure:"}︡{"stdout":" 15133\namended ff"}︡{"stdout":" -3\namended ff"}︡{"stdout":" 2*I - 1\nFound factor: 37 of 15133\nFound factor: 409 of 15133\n\nn: 15509\nFound factor: 1193 of 15509"}︡{"stdout":"\nFound factor: 13 of 15509\n\nn: 15529\namended ff"}︡{"stdout":" I + 2\nFound factor: 53 of 15529"}︡{"stdout":"\nFound factor: 293 of 15529\n\nn: 15613\nFound factor: 1201 of 15613"}︡{"stdout":"\nFound factor: 13 of 15613\n\nn: 15977\nFound factor: 13 of 15977"}︡{"stdout":"\n\nn: 16109\nFound factor: 181 of 16109"}︡{"stdout":"\nFound factor: 89 of 16109\n\nn: 16133\namended ff"}︡{"stdout":" I + 3\namended ff"}︡{"stdout":" I + 3\nFound factor: 949 of 16133\nFound factor: 221 of 16133\n\nn: 16601\nFound factor: 1277 of 16601"}︡{"stdout":"\nFound factor: 13 of 16601\n\nn: 16609\nFound factor: 17 of 16609"}︡{"stdout":"\nFound factor: 977 of 16609\n\nn: 16733\namended ff"}︡{"stdout":" 2*I + 1\namended ff"}︡{"stdout":" I - 2\nfailure:"}︡{"stdout":" 16733\nFound factor: 577 of 16733"}︡{"stdout":"\n\nn: 17177\namended ff"}︡{"stdout":" 2*I - 2\nFound factor: 193 of 17177\nFound factor: 89 of 17177\n\nn: 17261\namended ff"}︡{"stdout":" I - 3\nFound factor: 421 of 17261\nFound factor: 41 of 17261\n\nn: 17693\namended ff"}︡{"stdout":" 2*I - 1\nFound factor: 1361 of 17693\nFound factor: 13 of 17693\n\nn: 17797\namended ff"}︡{"stdout":" 3*I\nFound factor: 169 of 17797\nFound factor: 1369 of 17797\n\nn: 17849\nFound factor: 1373 of 17849"}︡{"stdout":"\nFound factor: 13 of 17849\n\nn: 17893\namended ff"}︡{"stdout":" 2*I + 1\namended ff"}︡{"stdout":" 2*I + 3\nFound factor: 617 of 17893\nFound factor: 29 of 17893\n\nn: 18317\namended ff"}︡{"stdout":" I + 3\nFound factor: 13 of 18317\nFound factor: 1409 of 18317\n\nn: 18409\nFound factor: 449 of 18409"}︡{"stdout":"\nFound factor: 41 of 18409\n\nn: 18989\nFound factor: 1117 of 18989"}︡{"stdout":"\nFound factor: 17 of 18989\n\nn: 19109\nFound factor: 97 of 19109"}︡{"stdout":"\nFound factor: 197 of 19109\n\nn: 19409\nFound factor: 13 of 19409"}︡{"stdout":"\n\nn: 19633\namended ff"}︡{"stdout":" 2*I - 2\nFound factor: 29 of 19633\nFound factor: 677 of 19633\n\nn: 19637\namended ff"}︡{"stdout":" 3*I + 1\nFound factor: 73 of 19637\nFound factor: 269 of 19637\n\nn: 20381\n"}︡{"stdout":"Found a factor in 38.95% of the cases\n"}︡{"stdout":"Failure List: [(3281, 'even'), (3281, 'odd')]\n"}︡{"stdout":"Elapsed time 120.31\n"}︡{"done":true}︡









