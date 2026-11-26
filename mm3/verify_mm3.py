"""verify_mm3.py — cálculo exacto (racional) del cociente r^2 = (distancia máx / distancia mín)^2
para n puntos en el espacio. Las coordenadas del archivo se toman como números
de punto flotante y se convierten exactamente a fracciones, así que el valor
impreso es el cociente verdadero de esa configuración (sin redondeo).

Uso: python3 verify_mm3.py archivo.txt [récord_publicado]
Formato: líneas "x y z" (las que empiezan por # se ignoran; si la primera línea
tiene 1-2 números se toma como encabezado).
"""
import sys
from fractions import Fraction
from decimal import Decimal, getcontext

getcontext().prec = 40


def main():
    pts = []
    with open(sys.argv[1]) as f:
        for line in f:
            p = line.split()
            if not p or line.startswith("#") or len(p) < 3:
                continue
            pts.append(tuple(Fraction(float(v)) for v in p[:3]))
    n = len(pts)
    dmin = dmax = None
    for i in range(n):
        for j in range(i + 1, n):
            d = sum((pts[i][k] - pts[j][k]) ** 2 for k in range(3))
            if d == 0:
                print("NO válido: puntos repetidos"); return 1
            dmin = d if dmin is None or d < dmin else dmin
            dmax = d if dmax is None or d > dmax else dmax
    r2 = dmax / dmin
    r2d = Decimal(r2.numerator) / Decimal(r2.denominator)
    print(f"n = {n}   r^2 = {r2d:.20f}   r = {r2d.sqrt():.20f}")
    if len(sys.argv) > 2:
        rec = sys.argv[2]
        R = Decimal(rec)
        ulp = Decimal(1).scaleb(-(len(rec) - rec.index(".") - 1))
        if r2d < R:
            print(f"MEJORA SEGURA: {r2d:.12f} < {rec} (el récord publicado es >= {rec})")
        elif r2d < R + ulp:
            print(f"DENTRO DEL REDONDEO: el récord está en [{rec}, {R + ulp}); hay que comparar con más cifras")
        else:
            print(f"No mejora: {r2d:.12f} >= {R + ulp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
