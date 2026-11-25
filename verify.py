"""verify.py — verificación EXACTA (aritmética racional) de un empaquetamiento
de n cuadrados unitarios en un triángulo equilátero.

Convierte cada cuadrado a coordenadas racionales exactas:
  - centro: racional (redondeado del archivo)
  - rotación: cos = (1-q^2)/(1+q^2), sin = 2q/(1+q^2) con q = tan(θ/2) racional,
    así cada cuadrado tiene lados de longitud EXACTAMENTE 1.
Luego comprueba, sin redondeos:
  - que ningún par de cuadrados se solape (teorema del eje separador),
  - que todos los vértices estén dentro del triángulo de lado s
    (las desigualdades con √3 se elevan al cuadrado, así que todo es racional).
Si la configuración del optimizador tiene solapamientos minúsculos (~1e-11),
la separa un poco expandiendo los centros y recalcula s.

Uso: python3 verify.py archivo.txt [record_actual] [dibujo.svg]
"""
import math
import sys
from fractions import Fraction as F


def load(path):
    with open(path) as f:
        n, s, _ = f.readline().split()
        n = int(n)
        sq = [tuple(map(float, f.readline().split())) for _ in range(n)]
    return n, float(s), sq


def rat_square(x, y, t, digits=14):
    q = F(math.tan(t / 2)).limit_denominator(10**digits)
    c = (1 - q * q) / (1 + q * q)
    s = 2 * q / (1 + q * q)
    cx = F(x).limit_denominator(10**digits)
    cy = F(y).limit_denominator(10**digits)
    h = F(1, 2)
    verts = []
    for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
        verts.append((cx + h * a * c - h * b * s, cy + h * a * s + h * b * c))
    return verts


def separated(P, Q):
    """True si los polígonos convexos P y Q tienen interiores disjuntos (exacto)."""
    for poly in (P, Q):
        m = len(poly)
        for k in range(m):
            x1, y1 = poly[k]
            x2, y2 = poly[(k + 1) % m]
            ax, ay = y2 - y1, x1 - x2  # normal al lado
            pa = [ax * x + ay * y for x, y in P]
            pb = [ax * x + ay * y for x, y in Q]
            if max(pa) <= min(pb) or max(pb) <= min(pa):
                return True
    return False


def inside(v, s):
    x, y = v
    # y >= 0 ; √3 x >= y ; √3 (s - x) >= y
    return y >= 0 and x >= 0 and 3 * x * x >= y * y and s - x >= 0 and 3 * (s - x) ** 2 >= y * y


def needed_side(polys):
    r3 = math.sqrt(3)
    pts = [v for p in polys for v in p]
    miny = min(float(y) for _, y in pts)
    minl = min(float(x) - float(y) / r3 for x, y in pts)
    maxr = max(float(x) + float(y) / r3 for x, y in pts)
    return maxr - minl - 2 * miny / r3, miny, minl


def build(sq, eps):
    cx = sum(p[0] for p in sq) / len(sq)
    cy = sum(p[1] for p in sq) / len(sq)
    return [((x - cx) * (1 + eps) + cx, (y - cy) * (1 + eps) + cy, t) for x, y, t in sq]


def check(sq):
    polys = [rat_square(*p) for p in sq]
    n = len(polys)
    for i in range(n):
        for j in range(i + 1, n):
            # descarte rápido por distancia
            if (sq[i][0] - sq[j][0]) ** 2 + (sq[i][1] - sq[j][1]) ** 2 >= 2.0:
                continue
            if not separated(polys[i], polys[j]):
                return None, (i, j)
    sf, miny, minl = needed_side(polys)
    # trasladar a la posición canónica con racionales y fijar s racional por encima
    r3 = math.sqrt(3)
    dy = -F(miny)
    dx = -F(minl + miny / r3)
    polys = [[(x + dx, y + dy) for x, y in p] for p in polys]
    for pad in (1e-13, 1e-12, 1e-11, 1e-10, 1e-9):
        s = F(sf + pad).limit_denominator(10**15)
        # pequeño margen extra hacia abajo/izquierda
        shift = F(pad / 4)
        P2 = [[(x + shift, y + shift / 4) for x, y in p] for p in polys]
        if all(inside(v, s) for p in P2 for v in p):
            return s, P2
    return None, "contención"


def svg(polys, s, path, title):
    sc = 500 / float(s)
    h = float(s) * math.sqrt(3) / 2
    W, H = 540, int(h * sc) + 60
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="20" y="22" font-family="sans-serif" font-size="15">{title}</text>']
    def P(x, y):
        return f"{20 + float(x) * sc:.2f},{H - 15 - float(y) * sc:.2f}"
    out.append(f'<polygon points="{P(0, 0)} {P(s, 0)} {P(float(s) / 2, h)}" fill="none" stroke="black" stroke-width="1.5"/>')
    for p in polys:
        out.append('<polygon points="' + " ".join(P(x, y) for x, y in p) + '" fill="#b9b9b9" stroke="black" stroke-width="0.8"/>')
    out.append("</svg>")
    with open(path, "w") as f:
        f.write("\n".join(out))


def main():
    path = sys.argv[1]
    record = float(sys.argv[2]) if len(sys.argv) > 2 else None
    n, s_claim, sq = load(path)
    for eps in [0, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]:
        s, res = check(build(sq, eps))
        if s is not None:
            sf = float(s)
            print(f"n={n}: VERIFICADO EXACTAMENTE con s = {sf:.10f} (expansión {eps:g})")
            if record is not None:
                d = record - sf
                print(f"  récord actual {record:.10f} -> {'MEJORA de %.3e' % d if d > 0 else 'no mejora (%.3e)' % d}")
            if len(sys.argv) > 3:
                svg(res, s, sys.argv[3], f"{n} cuadrados unitarios, s = {sf:.6f}")
                print(f"  dibujo: {sys.argv[3]}")
            return 0
    print(f"n={n}: NO se pudo verificar ({res})")
    return 1


if __name__ == "__main__":
    sys.exit(main())
