import math,random
_default_error = 1e-5
_max_iteration = 50
_default_decimals = 7
_default_width = 15

def print_defaults():
    print(f"error{_default_error}")
    print(f"max_iteration{_max_iteration}")


def _tabprint(l, decimals=None, width=None):
    if decimals is None:
        decimals = _default_decimals
    if width is None:
        width = _default_width
    def fmt(x):
        # Para float
        if isinstance(x, float):
            return f"{x:<{width}.{decimals}f}"
        
        # Para complex (viene de cmath, sympy.evalf(), etc.)
        if isinstance(x, complex):
            re = f"{x.real:.{decimals}f}"
            im = f"{abs(x.imag):.{decimals}f}"
            sign = "+" if x.imag >= 0 else "-"
            return f"{re}{sign}{im}i".ljust(width)
        
        # Para int, str, sympy.Expr, etc.
        return f"{str(x):{width}}"
    
    print("".join(fmt(x) for x in l))


def _default_iteration(fi, xinit, headers, cells, max_ite=_max_iteration,  eps=_default_error ):
    x_now = xinit
    iterations = 0
    rel_error = 100
    # table
    _tabprint(headers)
    while True:
        x_prev = x_now
        x_now = fi(x_now)
        iterations+=1
        _tabprint(cells(x_prev, x_now, iterations))
        if x_now != 0: 
            rel_error = abs( (x_now - x_prev) / x_now ) * 100
        if rel_error < eps or iterations >= max_ite:
            break
        #if (abs(x_now - x_prev) < eps or iterations >= max_ite): 
            #break
    return x_now, rel_error, iterations

#def binary_search(f, left, right, )

def bisection():
    print("No implemented")

def false_position():
    print("No implemented")


def newthon_rapshon(xinit, f, d1f, eps=_default_error, max_ite=_max_iteration   ):
    def cells(x_now, x_next, i):
        return [i, x_now, f(x_now), d1f(x_now), x_next ]

    return _default_iteration(fi=lambda x: x - f(x)/d1f(x), 
                       xinit=xinit, 
                       headers=["ite", "xi", "fx","d1f", "xi+1" ], 
                       cells=cells, 
                       eps=eps,
                       max_ite=max_ite)
    

def mod_newthon_raphson(xinit, f, d1f, d2f, eps=_default_error, max_ite=_max_iteration   ):
    def cells(x_now, x_next, i):
        return [i, x_now, f(x_now), d1f(x_now), d2f(x_now), x_next ]

    g = lambda x: x - (f(x) * d1f(x)) / (d1f(x)**2 - f(x)*d2f(x))
    return _default_iteration(fi=g, 
                       xinit=xinit, 
                       headers=["ite", "xi", "fx","d1f","d2f" ,"xi+1" ], 
                       cells=cells, 
                       eps=eps,
                       max_ite=max_ite)
    


def quadroot(r, s):
    disc = r**2 + 4*s

    t = math.sqrt(abs(disc))
    if disc > 0:
        r1 = (r + t) / 2
        r2 = (r - t) / 2
        i1 = i2 = 0 
    else:
        r1 = r2 = r / 2
        i1 = t / 2
        i2 = -i1
    return r1, i1, r2, i2

def bairstow(coef, r_init, s_init, eps=_default_error, max_ite=_max_iteration):
    # Copy
    coef = coef[:]
    r = r_init
    s = s_init
    iterations = 0
    r_error = 1
    s_error = 1
    degree = len(coef) - 1
    roots = [(0.0,0.0)] * (degree)
    status = 0
    b = [0.0] * (degree+1)
    c = [0.0] * (degree+1)

    while True:
        if degree < 3 or iterations >= max_ite: break
        iterations = 0
        print(f"Iteration with degree {degree}")
        while True:
            iterations+=1
            # division of polynomials
            b[degree] = coef[degree]
            b[degree-1] = coef[degree-1] + r * b[degree]
            c[degree] = b[degree]
            c[degree-1] = b[degree-1] + r * c[degree]
            # following coefficients 
            for i in range(degree-2, -1, -1):
                b[i] = coef[i] + r * b[i+1] + s * b[i+2]
                c[i] = b[i] + r * c[i+1] + s * c[i+2]
            det = c[2]**2 - c[3]*c[1]
            if det != 0:
                dr = (-b[1]*c[2]+b[0]*c[3]) / det
                ds = (-b[0]*c[2] + b[1]*c[1])/ det
                r+=dr
                s+=ds
                if r !=0: r_error = abs(dr / r)*100 
                else: r_error = float("inf")
                if s !=0: s_error = abs(ds / s)*100 
                else: s_error = float("inf")
            else:
                r+=1
                s+=1
                iterations = 0
            if (r_error <= eps and s_error <= eps) or iterations >= max_ite:
                break
        real1, imag1, real2, imag2 =  quadroot(r, s)
        #roots.append((real1, imag1))
        #roots.append((real2, imag2))
        roots[degree-1] = (real1, imag1)
        roots[degree-2] = (real2, imag2)
        degree -= 2
        # update coefficients
        for i in range(0, degree+1): coef[i] = b[i+2]
    if iterations < max_ite:
        if degree == 2: 
            r = -coef[1] / coef[2]
            s = -coef[0] / coef[2]
            real1, imag1, real2, imag2 =  quadroot(r, s)
            roots[degree-1] = (real1, imag1)
            roots[degree-2] = (real2, imag2)
        elif degree==1:
            roots[degree-1] = (-coef[0]/coef[1], 0)
    else:
        # there was an error
        print("does not converge")
        status=1
    return roots, status



if __name__ == "__main__":
    coe = [1, 3.5, 2.75, 2.125, 3.875, 125]

    def f(x): return x**5 + 3.5*x**4 + 2.75*x**3 + 2.125*x**2 + 3.875*x + 125

    r,s=bairstow(coe, -1, -1, eps=0.001)
    print()
    print()
    print()
    print()
    for i in r: print(i)

