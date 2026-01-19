import math
from math import exp, sqrt, erf, ceil
from math import factorial as fac

def gaussian_center_weight(sigma, t):
    """ Weight of the gaussian of std deviation s, on the interval [-t, t]
    :param x: (float)
    :param y: (float)
    :returns: erf( t / (sigma*\sqrt 2) )
    """
    return erf(t / (sigma * sqrt(2.)))


def binomial(x, y):
    """ Binomial coefficient
    :param x: (integer)
    :param y: (integer)
    :returns: y choose x
    """
    try:
        binom = fac(x) // fac(y) // fac(x - y)
    except ValueError:
        binom = 0
    return binom


def centered_binomial_pdf(k, x):
    """ Probability density function of the centered binomial law of param k at x
    :param k: (integer)
    :param x: (integer)
    :returns: p_k(x)
    """
    return binomial(2*k, x+k) / 2.**(2*k)


def build_centered_binomial_law(k):
    """ Construct the binomial law as a dictionnary
    :param k: (integer)
    :param x: (integer)
    :returns: A dictionnary {x:p_k(x) for x in {-k..k}}
    """
    D = {}
    for i in range(-k, k+1):
        D[i] = centered_binomial_pdf(k, i)
    return D

def build_centered_normal(sigma):
    """ Construct the normal law as a dictionary
    :param k: (integer)
    :param x: (integer)
    :returns: A dictionnary {x:p_k(x) for x in {-k..k}}
    """

    kmax = 1
    while True:
        partial_weight = 2*exp(-kmax**2/(2*sigma**2))
        if partial_weight < 2**(-400):
            break
        kmax += 1
    
    weight = 1
    for k in range(kmax, 0, -1):
        partial_weight = 2*exp(-k**2/(2*sigma**2))
        weight += partial_weight
        k += 1

    D = {}
    for i in range(-kmax, kmax+1):
        D[i] = exp(-i**2/(2*sigma**2))/weight

    return D


def mod_switch(x, q, rq):
    """ Modulus switching (rounding to a different discretization of the Torus)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    :param rq: output modulus (integer)
    """
    return int(round(1.* rq * x / q) % rq)


def mod_centered(x, q):
    """ reduction mod q, centered (ie represented in -q/2 .. q/2)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    """
    a = x % q
    if a < q/2:
        return a
    return a - q


def build_mod_switching_error_law(q, rq):
    """ Construct Error law: law of the difference introduced by switching from and back a uniform value mod q
    :param q: original modulus (integer)
    :param rq: intermediate modulus (integer)
    """
    D = {}
    V = {}
    for x in range(q):
        y = mod_switch(x, q, rq)
        z = mod_switch(y, rq, q)
        d = mod_centered(x - z, q)
        D[d] = D.get(d, 0) + 1./q
        V[y] = V.get(y, 0) + 1

    return D


def law_convolution(A, B):
    """ Construct the convolution of two laws (sum of independent variables from two input laws)
    :param A: first input law (dictionnary)
    :param B: second input law (dictionnary)
    """

    C = {}
    for a in A:
        for b in B:
            c = a+b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C


def law_product(A, B):
    """ Construct the law of the product of independent variables from two input laws
    :param A: first input law (dictionnary)
    :param B: second input law (dictionnary)
    """
    C = {}
    for a in A:
        for b in B:
            c = a*b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C

def clean_dist(A):
    """ Clean a distribution to accelerate further computation (drop element of the support with proba less than 2^-300)
    :param A: input law (dictionnary)
    """
    B = {}
    for (x, y) in A.items():
        if y>2**(-400):
            B[x] = y
    return B


def iter_law_convolution(A, i):
    """ compute the -ith forld convolution of a distribution (using double-and-add)
    :param A: first input law (dictionnary)
    :param i: (integer)
    """
    D = {0: 1.0}
    i_bin = bin(i)[2:]  # binary representation of n
    for ch in i_bin:
        print("iter", ch)
        D = law_convolution(D, D)
        D = clean_dist(D)
        if ch == '1':
            D = law_convolution(D, A)
            D = clean_dist(D)

    return D


def tail_probability(D, t):
    '''
    Probability that an drawn from D is strictly greater than t in absolute value
    :param D: Law (Dictionnary)
    :param t: tail parameter (integer)
    '''
    s = 0
    ma = max(D.keys())
    if t >= ma:
        return 0
    for i in reversed(range(int(ceil(t)), ma)):  # Summing in reverse for better numerical precision (assuming tails are decreasing)
        s += D.get(i, 0) + D.get(-i, 0)
    return s


def find_tail_for_probability(D, p):
    s = 0
    ma = max(D.keys())
    for i in reversed(range(1, ma)):  # Summing in reverse for better numerical precision (assuming tails are decreasing)
        s += D.get(i, 0) + D.get(-i, 0)
        if s >= p:
            return i-1
    return s

CBD3 = build_centered_binomial(3)  # for s and r, 2 for 256-bit security
CBD2 = build_centered_binomial(4)  # for e and e', 2 for 256-bit security
n, k = 256, 2           # k changed from 3 to 4 (256,5) for 256-bit security
q = 7681
d_v=6
d_u=12

S = CBD3
for _ in range(2):
    S = law_convolution(S, CBD3)
E_pub = CBD3
for _ in range(2):
    E_pub = law_convolution(E_pub, CBD3)
R = CBD3
for _ in range(1):
    R = law_convolution(R, CBD3)
Eprime = CBD2 # for u
for _ in range(1):
    Eprime = law_convolution(Eprime, CBD2)

E_extra = CBD2 # for v

E_v=build_mod_switching_error_law(q,2**d_v)

E_u=build_mod_switching_error_law(q,2**d_u)

# T1 = sum(ei) * r (n×k terms)
T1 = law_product(E_pub, R)
T2 = law_convolution(E_u, Eprime)
T3 = law_product(S, T2)
T4 = law_convolution(E_extra, E_v)
tot_terms_t1_t3 = n * k

# noise_subtotal = T1 - T3 (n×k terms)

# T1_sum=iter_law_convolution(T1, tot_terms_t1_t3)
# T2_sum=iter_law_convolution(T2, tot_terms_t1_t3)
# T3_sum=iter_law_convolution(T3, tot_terms_t1_t3)
sum_T1=iter_law_convolution(T1, tot_terms_t1_t3)
sum_T3=iter_law_convolution(T3, tot_terms_t1_t3)


# noise_subtotal = law_convolution(T1_sum, T2_sum)
noise_total = law_convolution(sum_T1, sum_T3)


# Finally combine: (T1 - T3) expanded over n×k + T2 expanded over n
noise_total = law_convolution(noise_total, T4)

# Target probability (consistent with your original script)
target_prob = 2**(-128) / n

# Calculate final tail bound
f_final = find_tail_for_probability(noise_total, target_prob)
print(f"Kyber-CPA tail-bound: {f_final}")

# Verify if modulus q satisfies q > 4 * f_final (check from your original script)

required_condition = q > 4 * f_final
print(f"is q={q} enough: {required_condition}")
if not required_condition:
    print(f"The requested q: {4 * f_final}")
    print(f"Recommendation: 12289, 18433, 40961")
else:
    print(f"Rest: {q/4 - f_final}")

print("log2 of k1*n * tail probability at (q/4):", math.log(n * tail_probability(noise_total, q/4), 2))

"""
128-bit security:
Kyber-CPA tail-bound: 1873
is q=7681 enough: True
Rest: 47.25
log2 of k1*n * tail probability at (q/4): -134.26362440826426

256-bit security
Kyber-CPA tail-bound: 1830
is q=7681 enough: True
Rest: 90.25
log2 of k1*n * tail probability at (q/4): -141.07202112243309
"""
