import numpy as np

PATH = "./assets/"
E = (2 ** 16) + 1


def fetch_messages():
    messages = []
    for ctr in range(1000):
        f = open(PATH + "msg_" + str(ctr) + ".txt", 'r')
        messages.append(int(f.read()))
        f.close()
    return messages


def fetch_traces():
    traces = np.zeros((1000, 40))
    for ctr in range(1000):
        f = open(PATH + "curve_" + str(ctr) + ".txt", 'r')
        trace = f.read().split(" ")
        trace.remove('')
        for k in range(40):
            traces[ctr][k] = float(trace[k])
        f.close()

    return traces


def fetch_modulo():
    f = open(PATH + "N.txt")
    N = f.read()
    f.close()
    return int(N)


def exponentiation(message, key_hyp, N):
    M = message
    for ctr in range(len(key_hyp) - 2, -1, -1):
        M = (M ** 2) % N
        if key_hyp[ctr] == 1:
            M = (M * message) % N
        else:
            if ctr == 0:
                M = (M ** 2) % N
    return bin(M).count('1')


def CPA(messages, traces, N):
    key_hyp = [1]
    hamming_weight_for_zero = np.zeros((1000, 1))
    hamming_weight_for_one = np.zeros((1000, 1))
    ctr = 1
    while ctr < len(traces[0]):
        for k in range(len(messages)):
            hamming_weight_for_zero[k] = exponentiation(messages[k], [0] + key_hyp, N)
            hamming_weight_for_one[k] = exponentiation(messages[k], [1] + key_hyp, N)
        if np.corrcoef(hamming_weight_for_zero, traces[:, ctr:ctr+1], False)[1][0] > \
                np.corrcoef(hamming_weight_for_one, traces[:, ctr:ctr+1], False)[1][0]:
            key_hyp = [0] + key_hyp
            ctr += 1
        else:
            key_hyp = [1] + key_hyp
            ctr += 2
    key_hyp.reverse()
    return key_hyp


def keyByFactorization(N):
    (p, q) = primeFactors(N)
    phi = (p - 1) * (q - 1)
    d = invmod(phi)
    return d


def primeFactors(N):
    i = 2
    while i**2 < N:
        if N % i:
            i += 1
        else:
            return i, N // i
    return False


def gcd(a, b):
    x, lastx, y, lasty = 0, 1, 1, 0
    while b != 0:
        a, (quotient, b) = b, divmod(a, b)
        x, lastx = lastx - quotient * x, x
        y, lasty = lasty - quotient * y, y
    return a, lastx * (-1 if a < 0 else 1), lasty * (-1 if b < 0 else 1)


def invmod(phi):
    g, x , y = gcd(E, phi)
    if g != 1:
        raise  ValueError
    return x % phi


N = fetch_modulo()
key_by_CPA = CPA(fetch_messages(), np.asarray(fetch_traces()), N)
key_string = ""
for i in key_by_CPA:
    key_string += str(i)

print("Key obtained by CPA : 0b" + key_string)

key_by_factorization = keyByFactorization(N)
key_string = bin(key_by_factorization)

print("Key obtained by N factorization :", key_string)
f = open("d_7.txt", 'w')
f.write(key_string)
f.close()