#!/usr/bin/python

R.<t> = PolynomialRing(QQ)
t = R.gens()[0]

memoize_m = {}
# number of set partitions that have a given shape
def m(pi):
    if tuple(pi) not in memoize_m:
        pi_t = list(Partition(pi).conjugate())
        pi_t.append(0)
        extra_factor = prod([factorial(pi_t[j] - pi_t[j+1]) for j in range(len(pi_t) - 1)])
        memoize_m[tuple(pi)] = multinomial(list(pi)) / extra_factor
    return memoize_m[tuple(pi)]

memoize_charpoly = {}
# pi is a *partition*, not a set partition
def CharacteristicPoly(pi):
    pi = tuple(pi)
    if pi not in memoize_charpoly:
        memoize_charpoly[pi] = prod([falling_factorial(t-1, part-1) for part in pi])
    return memoize_charpoly[pi]



# Kazhdan-Lusztig polynomial
memoize_KL = {1:1}
def KL(n):
    assert n >= 1
    if n not in memoize_KL:  # note: base case n=1 in dict already, so n >= 2
        result = 0
        for pi in Partitions(n):
            if len(pi) < n:
                result += m(pi) * CharacteristicPoly(pi) * KL(len(pi))
        # make sure it's a strictly lower degree than degree_bound
        degree_bound = floor(n/2)         
        memoize_KL[n] = -result.truncate(degree_bound) 
    return memoize_KL[n]

