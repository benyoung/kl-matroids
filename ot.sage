#!/usr/bin/python

base_ring.<q> = PolynomialRing(QQ)
Sym = SymmetricFunctions(base_ring)
p = Sym.powersum()
s = Sym.schur()
h = Sym.homogeneous()

def decompose_partition_into_rectangles(p):
    return [(i, list(p).count(i)) for i in reversed(range(max(p)+1)) if i in p]       

# M = (m1, m2, ... )
# iterate over all tuples of partitions of M
# this probably doesn't need to be its own function but whatever
def iterate_over_partition_tuples(M):
    return xmrange_iter([Partitions(u) for u in M])

def rep_dimension(f):
    d = f.monomial_coefficients()
    result = 0
    for mu in d:
        result += Partition(mu).dimension() * d[mu]
    return result

def Lie(n):
    result = Sym(0)
    n = Integer(n)
    for d in n.divisors():
        result += moebius(d) * p([d])**(n/d)
    return result/n

def multiplicity(mu, i):
    return mu.to_exp()[i-1]

def regular_representation(n):
    result = Sym(0)
    for llambda in Partitions(n):
        result += s(llambda) * StandardTableaux(llambda).cardinality()
    return result

def graded_dimension(rep):
    result = Sym(0)
    for (partition, coeff) in rep:
        result +=  StandardTableaux(partition).cardinality() * coeff
    return result

def B_coeff(n,j):
    result = 0
    for llambda in Partitions(n):
        if n-len(llambda) == j:
            prod = Sym(1)
            for i in range(1, max(llambda)+1):
                m_i = multiplicity(llambda, i)
                #print llambda, i, m_i
                factor = h([m_i]).plethysm(Lie(i))
                #print llambda, i, m_i, factor
                prod *= factor
            result += prod
    return result

def B(n):
    return sum(B_coeff(n,j) * q**j for j in range(n))

def R(n, degree_bound=9):
    result = 0
    for llambda in Partitions(n):
        coeff = base_ring(1)
        for (i, part) in enumerate(llambda):
            coeff *= q**(i * part)
        for (i,j) in llambda.cells():
            coeff *= (1-q**(degree_bound + j - i))
            hook = llambda.hook_length(i,j)
            coeff *= sum(q**(t*hook) for t in range(int(degree_bound/hook + 2)))
        coeff = ((1-q)*coeff).truncate(degree_bound)  
        result += coeff * s(llambda)
    return result

memoize_M = {
    1: s([1]),
    #2: s([2]),
    #3: s([3]) + q*s([1,1,1]),
    #4: s([4]) + q*s([2,1,1]) + q^2 * s([2,2]),
}
def M(n): # stub
    if n not in memoize_M:
        result = Sym(0)
        print n
        for i in range(n-1):
            result += q^i * M_coeff_from_OT(n, i)
        memoize_M[n] = result
    return memoize_M[n]

def extract_coeff(symm_func, i):
    result = Sym(0)
    for (part, old_poly) in symm_func:
        result += s(part) * old_poly[i] 
    return result

def M_compact_supp(n):
    result = Sym(0)
    Mn = M(n)
    for i in range(n):
        coeff = extract_coeff(Mn, i)
        result += coeff * q**(2*(n-1) - i)
    return result

#def OT(n, degree_bound=6):
#    return M(n).inner_tensor(R(n, degree_bound))

def M_coeff_from_OT(n, i):
    if n == 1:
        if i==0: 
            return s([1])
        else:
            return Sym(0)
    elif i > n-2:
        return Sym(0)
    else:
        result = extract_coeff(fake_OT(n), i)
        for k in range(1, i+1):
            left_piece = M_coeff_from_OT(n, i-k)
            right_piece = extract_coeff(R(n), k)
            result -= left_piece.inner_tensor(right_piece)
        return result

def fake_OT(n):
    result = Sym(0)
    for llambda in Partitions(n):
        if llambda != Partition([n]):
            L = decompose_partition_into_rectangles(llambda)
            M = [u[1] for u in L]
            B_rep = B(len(llambda))    
            for partition_list in iterate_over_partition_tuples(M):
                reference_rep = prod([s(mu) for mu in partition_list])
                term = reference_rep.scalar(B_rep)
                for i in range(len(partition_list)):
                    mu = partition_list[i]
                    r_i = L[i][0]
                    left_piece = M_compact_supp(r_i)
                    right_piece = R(r_i)
                    term *= s(mu).plethysm(left_piece.inner_tensor(right_piece))
                #    print llambda, partition_list, reference_rep.scalar(B_rep), reference_rep, left_piece.inner_tensor(right_piece), "->", term
                result += term

#            term = Sym(1)
#            for (r, m) in decompose_partition_into_rectangles(llambda):
#                term *= B(m).plethysm(R(r).inner_tensor(M_compact_supp(r)))
#            result += term
    return result

    
    


