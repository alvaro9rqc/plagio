from typing import Callable, Optional
import cmath

def muller_with_three(x0: complex, x1: complex, x2: complex,
                      eps: float, maxit: int, f: Callable[[complex], complex]):
    """
    Variante del método de Müller que recibe directamente x0, x1, x2.
    """
    iter_count = 0

    def fmt(z):
        if abs(z.imag) < 1e-12:
            return f"{z.real: .8f}"
        return f"({z.real: .6f}{z.imag:+.6f}j)"

    print(f"{'it':>3} | {'x0':>20} | {'x1':>20} | {'x2':>20} | {'xr(new)':>20} | {'dxr':>12} | {'|dxr|':>10} | {'Error %':>10} | {'f(x2)':>12}")
    print("-" * 150)

    xr_new = x2
    while True:
        iter_count += 1
        h0 = x1 - x0
        h1 = x2 - x1

        if h0 == 0 or h1 == 0:
            print("Error: h0 o h1 es cero (puntos coincidentes). Abortando.")
            break

        d0 = (f(x1) - f(x0)) / h0
        d1 = (f(x2) - f(x1)) / h1

        a = (d1 - d0) / (h1 + h0)
        b = a * h1 + d1
        c = f(x2)

        rad = cmath.sqrt(b * b - 4 * a * c)

        den = b + rad if abs(b + rad) > abs(b - rad) else b - rad
        if den == 0:
            print("Denominador cero al calcular dxr. Abortando.")
            break

        dxr = -2 * c / den
        xr_new = x2 + dxr

        # error relativo en porcentaje
        err_percent = (abs(dxr) / abs(xr_new) * 100) if xr_new != 0 else float("inf")

        print(f"{iter_count:3d} | {fmt(x0):>20} | {fmt(x1):>20} | {fmt(x2):>20} | {fmt(xr_new):>20} | {fmt(dxr):>12} | {abs(dxr):10.6e} | {err_percent:10.6f} | {fmt(c):>12}")

        #if abs(dxr) < eps * abs(xr_new) or iter_count >= maxit:
            #break
        if err_percent <= eps or iter_count >= maxit:
            break

        x0, x1, x2 = x1, x2, xr_new

    print("-" * 150)
    print(f"Iteraciones: {iter_count}, aproximación final: {fmt(xr_new)}")
    return xr_new


def muller(xr: Optional[complex] = None, h: Optional[float] = None,
           x0: Optional[complex] = None, x1: Optional[complex] = None, x2: Optional[complex] = None,
           eps: float = 1e-6, maxit: int = 50,
           f: Callable[[complex], complex] = lambda x: x):
    """
    Wrapper flexible:
    - Si se pasa xr y h => usa pseudocódigo (xr central, h escala).
    - Si se pasan x0, x1, x2 => usa esos valores iniciales.
    """
    if xr is not None and h is not None:
        x2 = complex(xr)
        x1 = complex(xr + h * xr)
        x0 = complex(xr - h * xr)
        return muller_with_three(x0, x1, x2, eps, maxit, f)
    elif x0 is not None and x1 is not None and x2 is not None:
        return muller_with_three(x0, x1, x2, eps, maxit, f)
    else:
        raise ValueError("Debes dar (xr y h) o bien (x0, x1, x2)")

"""

def f(x): return x**3 - 13*x-12

import muller

muller.muller_with_three(4.5, 5.5, 5, 1e-5, 50, f)

"""
