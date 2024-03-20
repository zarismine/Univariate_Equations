# Majid Mumtaz'method
from Crypto.Util.number import *
from sage.matrix.matrix2 import Matrix
def resultant(f1, f2, var):
    return Matrix.determinant(f1.sylvester_matrix(f2, var))

def coron(pol, X, Y, Z, l = 2, debug=False):
    if pol.nvariables() != 3:
        raise ValueError("pol is not bivariate")

    P.<x,y,z> = PolynomialRing(ZZ)
    pol = pol(x,y,z)
    xoffset = 0
    while pol(xoffset,0,0) == 0:
        xoffset += 1
    pol = pol(x+xoffset,y,z)
    while gcd(pol(0,0,0), X) != 1:
        X = next_prime(X, proof=False)

    while gcd(pol(0,0,0), Y) != 1:
        Y = next_prime(Y, proof=False)
    while gcd(pol(0,0,0), Z) != 1:
        Z = next_prime(Z, proof=False)

    pol = P(pol/gcd(pol.coefficients())) # seems to be helpful
    p00 = pol(0,0,0)
    delta = max(pol.degree(x),pol.degree(y),pol.degree(z)) # maximum degree of any variable

    W = max(abs(i) for i in pol(x*X,y*Y,z*Z).coefficients())
    u = W + ((1-W) % abs(p00))
    N = u*(X*Y*Z)^l # modulus for polynomials

    # Construct polynomials
    p00inv = inverse_mod(p00,N)
    polq = P(sum((i*p00inv % N)*j for i,j in zip(pol.coefficients(),pol.monomials())))
#     print(polq)
    polynomials = []
    for i in range(delta+l+1):
        for j in range(delta+l+1):
            for k in range(delta+l+1):
                if 0 <= i <= l and 0 <= j <= l and 0 <= k <= l:
                    polynomials.append(polq * x^i * y^j *z^k * X^(l-i) * Y^(l-j) * Z^(l-k))
                else:
                    polynomials.append(x^i * y^j * z^k * N)

    # Make list of monomials for matrix indices
    monomials = []
    for i in polynomials:
        for j in i.monomials():
            if j not in monomials:
                monomials.append(j)
    monomials.sort()
    L = matrix(ZZ,len(monomials))
    for i in range(len(monomials)):
        for j in range(len(monomials)):
            L[i,j] = polynomials[i](X*x,Y*y,Z*z).monomial_coefficient(monomials[j])

    # makes lattice upper triangular
    # probably not needed, but it makes debug output pretty
    L = matrix(ZZ,sorted(L,reverse=True))
    L = L.LLL()
    for i in range(64):
        pol1 = P(sum(map(mul, zip(L[i],monomials)))(x/X,y/Y,z/Z))
        if (gcd(pol1,pol))==1 and pol1(d,kc1,kc2)==0:
            print(i)
    pol1 = P(sum(map(mul, zip(L[0],monomials)))(x/X,y/Y,z/Z))
    pol2 = P(sum(map(mul, zip(L[1],monomials)))(x/X,y/Y,z/Z))
#     print(pol(d,kc1,kc2),pol2(d,kc1,kc2),pol3(d,kc1,kc2),pol4(d,kc1,kc2),pol5(d,kc1,kc2))
    h1 = resultant(pol, pol1, z)
    h2 = resultant(pol, pol2, z) #b,c
    ff = resultant(h1, h2, y)
    
    roots = ff.univariate_polynomial().roots()
    print(roots)



p =  4281528017348910901488388149849326033527889641499112395113403771496406064058889328362117032041400534222317207495790524164270193113426166773406205554981
q =  3322901489097492093478695211438889459420740749118370722694272573557370858396069565784465520956175806272428909647017798466127814899951252410631104587561
a =  1328109911095185674159095350341162669564591890
b =  1030748457882552013053809028358749402940050580
g =  1611887683986289320654889897425606366477748495956317632048294707249638258704676438393846077996717839759041


n = p*q

t = 130
d = getPrime(t)

e = ZZ(inverse(d, 2 * g * a * b))
ee = t+500-350
X, Y, Z = 2^t,2^ee,2^ee
PR.<x,y,z> = PolynomialRing(ZZ)
f = e^2*x^2+e*x*(y+z-2)-(y+z-1)-(n-1)*y*z
# print(f)
kc1 = (e*d - 1)//(q-1)
assert (e*d - 1) == kc1*(q-1)
kc2 = (e*d - 1)//(p-1)
# assert (e*d - 1) == kc2*(p-1)
print(len(bin(kc1))-2)
coron(f, X, Y, Z,l = 1)
