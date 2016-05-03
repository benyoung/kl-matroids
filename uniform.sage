#!/usr/bin/python

R.<t> = PolynomialRing(QQ)

# Kazhdan-Lusztig polynomial
memoize_KL = {(1,1):R(1)}
def KL(m,d):
    assert m >= 1
    assert d >= 1
    if (m,d) not in memoize_KL:  # note: base case m=1, d=1 in dict already, so m >= 2
        result=0
        for i in range(d):
            result += ((-1)^(i))*(binomial(m+d,i))*((t**(d-i))-1)
        for s in range(1,d):
            result += binomial(m+d,s) * ((t-1)**s) * KL(m,d-s)
                # make sure it's a strictly lower degree than degree_bound
        degree_bound = floor((d-1)/2)+1         
        memoize_KL[(m,d)] = -result.truncate(degree_bound) 
    return memoize_KL[(m,d)]

