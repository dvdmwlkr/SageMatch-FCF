︠e8f8895c-ddc9-4bcb-85b3-6b59f201b731s︠
import time

vers = version()
print vers; print

def PrimeDivisorCount(n):
    pd = prime_divisors(n)
    return len(pd)

#F = (x for x in range(15,5001000,2) if (PrimeDivisorCount(x) > 1) and not is_square(x) and not is_prime(x) and (x %3 != 0) and (x%5 != 0))
F = (x for x in range(17,20050,2) if not is_square(x) and not is_prime(x) and (x %3 != 0) and (x%5 != 0))
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
            print "Even Failures: \n  {0}   \nOdd Failures: \n  {1}".format(FailuresEvenLen, FailuresOddLen)
        ix2 += 1
    last_n = n
︡f2efaaef-a2ec-4f1d-8798-e1128aca70a7︡{"stdout":"SageMath version 7.3, Release Date: 2016-08-04\n\n"}︡{"stdout":"N:    1001   count:     93   failuresEvenLen:     9(9.68%)   failuresOddLen:     2(2.15%)   Success: 88.17%   time: 0.07\nN:    2009   count:    223   failuresEvenLen:    18(8.07%)   failuresOddLen:     4(1.79%)   Success: 90.13%   time: 0.14"}︡{"stdout":"\nN:    5017   count:    651   failuresEvenLen:    47(7.22%)   failuresOddLen:    15(2.30%)   Success: 90.48%   time: 0.39"}︡{"stdout":"\nN:   10001   count:   1415   failuresEvenLen:    90(6.36%)   failuresOddLen:    28(1.98%)   Success: 91.66%   time: 0.96"}︡{"stdout":"\nN:   20003   count:   3038   failuresEvenLen:   173(5.69%)   failuresOddLen:    70(2.30%)   Success: 92.00%   time: 2.32"}︡{"stdout":"\nEven Failures: \n  [119, 287, 343, 527, 623, 731, 779, 803, 959, 1003, 1127, 1331, 1343, 1411, 1519, 1679, 1687, 1691, 2023, 2123, 2147, 2263, 2599, 2603, 2651, 2747, 2807, 2863, 3239, 3247, 3403, 3479, 3587, 3647, 3667, 3707, 3791, 3827, 3971, 4151, 4223, 4247, 4319, 4487, 4499, 4579, 4763, 5027, 5327, 5339, 5371, 5383, 5543, 5627, 6167, 6191, 6239, 6319, 6347, 6443, 6499, 6503, 6559, 6611, 6727, 6859, 6887, 7223, 7339, 7343, 7403, 7567, 7571, 7619, 8023, 8227, 8279, 8459, 8483, 8519, 8531, 8651, 8687, 8899, 8927, 9023, 9079, 9247, 9407, 9799, 10199, 10291, 10307, 10327, 10423, 10447, 10727, 10747, 10823, 10951, 10999, 11023, 11099, 11207, 11263, 11303, 11327, 11543, 11711, 11723, 11879, 12167, 12179, 12319, 12431, 12587, 12767, 13067, 13123, 13223, 13379, 13391, 13439, 13739, 13823, 13943, 14003, 14191, 14267, 14507, 14603, 14611, 14903, 14959, 15127, 15479, 15623, 15967, 16019, 16199, 16379, 16439, 16571, 16639, 16643, 16847, 17063, 17083, 17111, 17311, 17363, 17639, 17687, 17867, 18071, 18227, 18599, 18631, 18767, 18779, 18983, 18991, 19103, 19307, 19439, 19547, 19567, 19723, 19771, 19871, 19879, 19883, 19967]   \nOdd Failures: \n  [533, 697, 1313, 1853, 2197, 2249, 2813, 3029, 3281, 3293, 3341, 3653, 3973, 4469, 4913, 5213, 5353, 5629, 5713, 5837, 5933, 6253, 6893, 8321, 8653, 9197, 9773, 9953, 10397, 11029, 11141, 11401, 11453, 11509, 11629, 11773, 11849, 12053, 12389, 12461, 12557, 12629, 12773, 13549, 13801, 14089, 14453, 15133, 15509, 15529, 15613, 15977, 16109, 16133, 16601, 16609, 16733, 17177, 17261, 17693, 17797, 17849, 17893, 18317, 18409, 18989, 19109, 19409, 19633, 19637]\n"}︡{"done":true}︡









