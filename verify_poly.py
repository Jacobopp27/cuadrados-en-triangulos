"""verify_poly.py — verificación rigurosa de un empaquetamiento de n cuadrados
unitarios dentro de un polígono regular de K lados (archivo de polypack).

- Cada cuadrado se convierte a coordenadas RACIONALES exactas, con lados de
  longitud exactamente 1 (rotación racional: cos=(1-q²)/(1+q²), sin=2q/(1+q²)).
- Que no haya superposiciones se comprueba con aritmética racional EXACTA
  (teorema del eje separador).
- La contención involucra cos/sen de múltiplos de 360°/K (irracionales): se
  calcula con 80 dígitos decimales (error < 1e-70) y el lado se redondea hacia
  ARRIBA con un margen de 1e-30, muy por encima de cualquier error numérico.

Uso: python3 verify_poly.py archivo.txt [record_actual] [dibujo.svg]
"""
import math
import sys
from decimal import Decimal as D, getcontext
from fractions import Fraction as F

getcontext().prec = 80


def dpi():
    # Machin: pi = 16 atan(1/5) - 4 atan(1/239)
    def atan_inv(x):
        x = D(x); x2 = x * x; term = 1 / x; s = term; k = 1; sign = -1
        while True:
            term /= x2
            t = term / (2 * k + 1)
            if t < D(10) ** -85:
                return s
            s += sign * t; sign = -sign; k += 1
    return 16 * atan_inv(5) - 4 * atan_inv(239)


PI = dpi()


def dcos(x):
    s = D(1); term = D(1); k = 0
    while True:
        k += 2
        term *= -x * x / (k * (k - 1))
        if abs(term) < D(10) ** -85:
            return s
        s += term


def dsin(x):
    s = x; term = x; k = 1
    while True:
        k += 2
        term *= -x * x / (k * (k - 1))
        if abs(term) < D(10) ** -85:
            return s
        s += term


def load(path):
    with open(path) as f:
        head = f.readline().split()
        n, s = int(head[0]), float(head[1])
        K = int(head[3]) if len(head) > 3 else 3
        sq = [tuple(map(float, f.readline().split())) for _ in range(n)]
    return n, K, s, sq


def rat_square(x, y, t, digits=14):
    q = F(math.tan(t / 2)).limit_denominator(10**digits)
    c = (1 - q * q) / (1 + q * q)
    s = 2 * q / (1 + q * q)
    cx = F(x).limit_denominator(10**digits)
    cy = F(y).limit_denominator(10**digits)
    h = F(1, 2)
    return [(cx + h * a * c - h * b * s, cy + h * a * s + h * b * c) for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))]


def separated(P, Q):
    for poly in (P, Q):
        m = len(poly)
        for k in range(m):
            x1, y1 = poly[k]
            x2, y2 = poly[(k + 1) % m]
            ax, ay = y2 - y1, x1 - x2
            pa = [ax * x + ay * y for x, y in P]
            pb = [ax * x + ay * y for x, y in Q]
            if max(pa) <= min(pb) or max(pb) <= min(pa):
                return True
    return False


def todec(fr):
    return D(fr.numerator) / D(fr.denominator)


def needed_side(polys, K):
    """Lado mínimo del K-gono regular (lado inferior horizontal, centro libre)
    que contiene todos los vértices: programa lineal resuelto enumerando
    tripletas de lados activos, en 80 dígitos."""
    normals = []
    for k in range(K):
        phi = -PI / 2 + 2 * PI * k / K
        normals.append((dcos(phi), dsin(phi)))
    pts = [(todec(x), todec(y)) for p in polys for x, y in p]
    h = [max(nx * x + ny * y for x, y in pts) for nx, ny in normals]
    best = None
    for p in range(K):
        for q in range(p + 1, K):
            for r in range(q + 1, K):
                ids = (p, q, r)
                M = [[normals[i][0], normals[i][1], D(1)] for i in ids]
                v = [h[i] for i in ids]

                def det3(A):
                    return (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
                            - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                            + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))
                d = det3(M)
                if abs(d) < D(10) ** -30:
                    continue
                sol = []
                for col in range(3):
                    Mc = [row[:] for row in M]
                    for t in range(3):
                        Mc[t][col] = v[t]
                    sol.append(det3(Mc) / d)
                cx, cy, a = sol
                if all(h[k] - normals[k][0] * cx - normals[k][1] * cy <= a + D(10) ** -60 for k in range(K)):
                    if best is None or a < best[0]:
                        best = (a, cx, cy)
    a, cx, cy = best
    # el centro (cx, cy) con apotema a contiene todo; holgura extra de 1e-30
    a += D(10) ** -30
    tan_pk = dsin(PI / K) / dcos(PI / K)
    return 2 * a * tan_pk, cx, cy


def build(sq, eps):
    cx = sum(p[0] for p in sq) / len(sq)
    cy = sum(p[1] for p in sq) / len(sq)
    return [((x - cx) * (1 + eps) + cx, (y - cy) * (1 + eps) + cy, t) for x, y, t in sq]


def check(sq, K):
    polys = [rat_square(*p) for p in sq]
    n = len(polys)
    for i in range(n):
        for j in range(i + 1, n):
            if (sq[i][0] - sq[j][0]) ** 2 + (sq[i][1] - sq[j][1]) ** 2 >= 2.0:
                continue
            if not separated(polys[i], polys[j]):
                return None, (i, j)
    s, cx, cy = needed_side(polys, K)
    return (s, cx, cy), polys


def svg(polys, K, s, cx, cy, path, title):
    s = float(s); cx = float(cx); cy = float(cy)
    R = s / (2 * math.sin(math.pi / K))
    verts = []
    for k in range(K):
        # vértices entre normales consecutivas
        ang = -math.pi / 2 + 2 * math.pi * (k + 0.5) / K
        verts.append((cx + R * math.cos(ang), cy + R * math.sin(ang)))
    xs = [v[0] for v in verts]; ys = [v[1] for v in verts]
    sc = 480 / max(max(xs) - min(xs), max(ys) - min(ys))
    W, H = 540, int((max(ys) - min(ys)) * sc) + 70

    def P(x, y):
        return f"{30 + (float(x) - min(xs)) * sc:.2f},{H - 15 - (float(y) - min(ys)) * sc:.2f}"
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="20" y="22" font-family="sans-serif" font-size="15">{title}</text>',
           '<polygon points="' + " ".join(P(x, y) for x, y in verts) + '" fill="none" stroke="black" stroke-width="1.5"/>']
    for p in polys:
        out.append('<polygon points="' + " ".join(P(x, y) for x, y in p) + '" fill="#b9b9b9" stroke="black" stroke-width="0.8"/>')
    out.append("</svg>")
    with open(path, "w") as f:
        f.write("\n".join(out))


def main():
    path = sys.argv[1]
    record = D(sys.argv[2]) if len(sys.argv) > 2 else None
    n, K, s_claim, sq = load(path)
    res = None
    for eps in [0, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]:
        r, res = check(build(sq, eps), K)
        if r is not None:
            s, cx, cy = r
            print(f"K={K} n={n}: VERIFICADO con s = {s:.12f} (expansión {eps:g})")
            if record is not None:
                d = record - s
                print(f"  récord actual {record} -> {'MEJORA de %.3e' % d if d > 0 else 'no mejora (%.3e)' % d}")
            if len(sys.argv) > 3:
                svg(res, K, s, cx, cy, sys.argv[3], f"{n} cuadrados en {K}-gono, s = {float(s):.6f}")
                print(f"  dibujo: {sys.argv[3]}")
            return 0
    print(f"K={K} n={n}: NO se pudo verificar (solapamiento {res})")
    return 1


if __name__ == "__main__":
    sys.exit(main())
