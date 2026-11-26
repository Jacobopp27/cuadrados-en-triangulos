"""digitize_poly.py — lee una imagen de récord del Packing Center (cuadrados
rellenos dentro de un polígono regular de K lados con el lado inferior
horizontal) y estima centro y ángulo de cada cuadrado, en el sistema de
polypack (polígono centrado en el origen, y hacia arriba).

Uso: python3 digitize_poly.py imagen K s_record n_esperado salida.txt
"""
import math
import sys
from collections import deque

import numpy as np
from PIL import Image


def components(mask):
    H, W = mask.shape
    seen = np.zeros_like(mask)
    comps = []
    for y0, x0 in zip(*np.nonzero(mask)):
        if seen[y0, x0]:
            continue
        q = deque([(y0, x0)]); seen[y0, x0] = True; pix = []
        while q:
            y, x = q.popleft(); pix.append((x, y))
            for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                yy, xx = y + dy, x + dx
                if 0 <= yy < H and 0 <= xx < W and mask[yy, xx] and not seen[yy, xx]:
                    seen[yy, xx] = True; q.append((yy, xx))
        comps.append(np.array(pix, dtype=float))
    return comps


def main():
    path, K, s_rec, n_exp, out = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    im = np.asarray(Image.open(path).convert("RGB")).astype(int)
    lum = im.mean(axis=2)
    # color de relleno = color no blanco más frecuente
    cols, cnts = np.unique(im.reshape(-1, 3), axis=0, return_counts=True)
    order = np.argsort(-cnts)
    fc = next(cols[i] for i in order if cols[i].mean() < 250 and cols[i].mean() > 120)
    fill = np.abs(im - fc).max(axis=2) <= 12
    # bordes: claramente más oscuros que el relleno
    dark = lum < min(fc.mean() - 45, 170)

    # contorno: proyección extrema de los píxeles oscuros sobre cada normal
    # (en la imagen y crece hacia abajo)
    ys, xs = np.nonzero(dark)
    normals = [(math.cos(-math.pi / 2 + 2 * math.pi * k / K), -math.sin(-math.pi / 2 + 2 * math.pi * k / K)) for k in range(K)]
    h = [max(nx * xs + ny * ys) for nx, ny in normals]
    # mínimos cuadrados: n_k·c + a = h_k
    A = np.array([[nx, ny, 1.0] for nx, ny in normals]); sol = np.linalg.lstsq(A, np.array(h), rcond=None)[0]
    cx, cy, a_px = sol
    # grosor del borde: se resta medio grosor aproximado (1 px)
    a_px -= 1.0
    CA = 1 / (2 * math.tan(math.pi / K))
    scale = a_px / (CA * s_rec)

    squares = []
    unit = scale * scale
    for P in components(fill):
        k = len(P) / unit
        if k < 0.35:
            continue
        cnt = max(1, round(k / 0.85))
        best = None
        for a10 in range(0, 900, 5):
            a = math.radians(a10 / 10)
            u = P[:, 0] * math.cos(a) + P[:, 1] * math.sin(a)
            v = -P[:, 0] * math.sin(a) + P[:, 1] * math.cos(a)
            area = (u.max() - u.min()) * (v.max() - v.min())
            if best is None or area < best[0]:
                best = (area, a, (u.max() + u.min()) / 2, (v.max() + v.min()) / 2, u.max() - u.min(), v.max() - v.min())
        _, a, cu, cv, lu, lv = best
        px = cu * math.cos(a) - cv * math.sin(a)
        py = cu * math.sin(a) + cv * math.cos(a)
        if cnt == 1:
            centers = [(px, py)]
        else:
            L = max(lu, lv); ax = (math.cos(a), math.sin(a)) if lu >= lv else (-math.sin(a), math.cos(a))
            centers = [(px + ax[0] * (-L / 2 + L * (i + 0.5) / cnt), py + ax[1] * (-L / 2 + L * (i + 0.5) / cnt)) for i in range(cnt)]
        t = (-a) % (math.pi / 2)
        for qx, qy in centers:
            squares.append(((qx + 0.5 - cx) / scale, (cy - qy - 0.5) / scale, t))
    print(f"{path}: {len(squares)} cuadrados (esperados {n_exp}), escala {scale:.2f} px/unidad")
    with open(out, "w") as f:
        f.write(f"{len(squares)} {s_rec:.6f} 0 {K}\n")
        for x, y, t in squares:
            f.write(f"{x:.6f} {y:.6f} {t:.6f}\n")
    return 0 if len(squares) == n_exp else 2


if __name__ == "__main__":
    sys.exit(main())
