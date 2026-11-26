"""verify_isos.py — verificación exacta (enteros) de un conjunto de puntos de la
cuadrícula n×n sin triángulos isósceles (problema 6.59 de AlphaEvolve).

Condición: para tres puntos distintos a, b, c del conjunto nunca se cumple
d(a,b) = d(a,c). Equivale a: desde cada punto, las distancias al cuadrado hacia
todos los demás son distintas entre sí (esto incluye el caso degenerado del
punto medio de un segmento).

Uso: python3 verify_isos.py n archivo.txt   (una línea "x y" por punto)
"""
import sys


def main():
    n = int(sys.argv[1])
    pts = []
    with open(sys.argv[2]) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                x, y = map(int, line.split()[:2])
                pts.append((x, y))
    assert len(set(pts)) == len(pts), "puntos repetidos"
    assert all(0 <= x < n and 0 <= y < n for x, y in pts), "punto fuera de la cuadrícula"
    for i, (ax, ay) in enumerate(pts):
        seen = {}
        for j, (bx, by) in enumerate(pts):
            if i == j:
                continue
            d = (ax - bx) ** 2 + (ay - by) ** 2
            if d in seen:
                print(f"NO válido: {pts[i]} equidista de {pts[seen[d]]} y {pts[j]}")
                return 1
            seen[d] = j
    print(f"VÁLIDO: {len(pts)} puntos en la cuadrícula {n}x{n}, sin triángulos isósceles")
    return 0


if __name__ == "__main__":
    sys.exit(main())
