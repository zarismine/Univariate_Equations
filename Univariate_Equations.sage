from sage.matrix.matrix2 import Matrix
from itertools import product
from re import findall
from subprocess import check_output
def flatter(M):
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    ret = check_output(["flatter"], input=z.encode())
    return matrix(M.nrows(), M.ncols(), map(int, findall(b"-?\\d+", ret)))
def resultant(f1, f2, var):
    return Matrix.determinant(f1.sylvester_matrix(f2, var))
def polynomials_con(F,N,F_exp,N_exp):
    P = N^N_exp
    for i,j in zip(F,F_exp):
        P *= i^j
    return P
def Univariate_Equations(N,F,r,R,Gamma,Eta,t):
    """
    params:
    N,r : N = 0 mod p^r
    F : Univariate polynomials set
    R : set of R[i] where F[i] = 0 mod p^R[i]
    Gamma : x[i] < N^gamma[i]
    Eta : p > N^eta
    t : special param about matrix-lattice dim
    return : x[i]
    """
    for i in range(len(F)):
        if Gamma[i]/R[i] > Eta:
            F.remove(F[i])
            Gamma.remove(Gamma[i])
            R.remove(R[i])
    n = len(F)
    S = [i for i in range(20)]
    assert len(F) == len(Gamma) == len(R)
    if prod(Gamma)/(r*prod(R)) > Eta^(n+1) or Eta < 1/sqrt(log(N,2)):
        return [0] * n
    bounds = [int(ceil(N^i)) for i in Gamma]
    G = Sequence([], F[0].parent())
    # for more equations,change this
    for item in product(S,repeat = n):
        if sum([i*j for i,j in zip(item,Gamma)]) < t*Eta:
            temp = t - sum([i*j for i,j in zip(R,item)])
            N_exp = max(ceil(temp),0)
            G.append(polynomials_con(F,N,item,N_exp))
    B, monomials = G.coefficient_matrix()
    monomials = vector(monomials)
    # print(monomials,factors)
    factors = [monomial(*bounds) for monomial in monomials]
    for i, factor in enumerate(factors):
        B.rescale_col(i, factor)
    print("lattice dim:",B.ncols(),B.nrows())
    print("start LLL")
    # try to flatter
#     B = B.dense_matrix().LLL()
    B = flatter(B)
    print("end LLL")
    B = B.change_ring(QQ)
    for i, factor in enumerate(factors):
        B.rescale_col(i, 1 / factor)
    R = B * monomials
    # change this if have more vars
    h1 = R[0]
    h2 = R[1]
    h3 = resultant(h1,h2,y)
    rx = h3.univariate_polynomial().roots()[0][0]
    ry = h1(rx,y).univariate_polynomial().roots()[0][0]
    return (rx,ry)

N = 0xbe9ccc83003bedf45421b58377b946f87dfd85be82124dc5d732070d77ef68e0231c3f34dc803a8984de0573db6d83ccea0bd53a885059a10cfa3764c658c4d42c5fa90ecad8573fff8f2c41e513278c59121e42ad83310fb22b4d20e7ada42c76f08891f38c92a1b1aac712bfa7d717a4c4802ed023f12c768972ca1b
e = 0x5dc97ed7250e57ce6fac4f57885c0538b1ea540fbaca79730470b6b990f7e861adc4c5fee3acdcd9ae9a2834b606ddfae01ade33edfa96a47a0ffc0036a4497a84c38b7cdac20c38f

E = ZZ(inverse_mod(e, N-1))
PR.<x,y> = PolynomialRing(ZZ)
f1 = E - x
f2 = N - y
Univariate_Equations(N-1,[f1,f2],1,[1,2],[0.25,0.5],0.42,10)