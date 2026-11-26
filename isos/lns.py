"""lns.py — reconstrucción EXACTA de una región de una solución sin isósceles.

Se fija todo lo que está fuera de la región, y un solucionador de restricciones
(CP-SAT de OR-Tools) busca la mayor cantidad de órbitas (simetría D2: grupos de 4
puntos reflejados en ambos ejes) dentro de la región sin crear ningún isósceles.
Si demuestra el óptimo, sabemos con certeza cuántos puntos caben ahí.

La región se describe en el cuadrante inferior izquierdo (0 <= x, y < n/2) con
cajas "x0:x1,y0:y1" (extremo derecho excluido); se aplica a las 4 esquinas por simetría.

Uso: python lns.py n solucion.txt "x0:x1,y0:y1[;...]" segundos [salida.txt]
"""
import sys, time
from collections import defaultdict
import numpy as np
from ortools.sat.python import cp_model


def load(fn):
    return [tuple(map(int, l.split()[:2])) for l in open(fn) if l.strip() and not l.startswith('#')]


def main():
    n = int(sys.argv[1]); P = set(load(sys.argv[2]))
    boxes = [tuple(map(int, b.replace(',', ':').split(':'))) for b in sys.argv[3].split(';')]
    secs = float(sys.argv[4]); out = sys.argv[5] if len(sys.argv) > 5 else None
    assert n % 2 == 0, "este modelo supone n par (órbitas de 4 casillas)"
    h = n // 2
    rep = lambda x, y: (min(x, n - 1 - x), min(y, n - 1 - y))
    cells = lambda r: [(x, y) for x in (r[0], n - 1 - r[0]) for y in (r[1], n - 1 - r[1])]

    R = {(x, y) for (x0, x1, y0, y1) in boxes for x in range(x0, min(x1, h)) for y in range(y0, min(y1, h))}
    F = sorted(p for p in P if rep(*p) not in R)            # puntos fijos
    old = len({rep(*p) for p in P if rep(*p) in R})         # órbitas actuales en la región
    Fa = np.array(F); m = len(F)
    Dmax = 2 * (n - 1) ** 2 + 1

    # T[i, d] = True si el punto fijo i ya tiene otro punto fijo a distancia^2 d
    T = np.zeros((m, Dmax), dtype=bool)
    for i in range(m):
        d = ((Fa - Fa[i]) ** 2).sum(1); d[i] = 0
        T[i, d[d > 0]] = True

    # 1) descartar órbitas de la región que ya chocan con los puntos fijos
    t0 = time.time(); V = []
    for r in sorted(R):
        K = cells(r); ok = True
        dK = [((Fa - np.array(k)) ** 2).sum(1) for k in K]
        for a, k in enumerate(K):
            others = list(dK[a]) + [(k[0] - q[0]) ** 2 + (k[1] - q[1]) ** 2 for b, q in enumerate(K) if b != a]
            if len(set(others)) < len(others): ok = False; break              # k como vértice
            if T[np.arange(m), dK[a]].any(): ok = False; break                # base (k, fijo) con vértice fijo
            for b in range(a + 1, len(K)):
                if (dK[a] == dK[b]).any(): ok = False; break                  # base (k, k') con vértice fijo
            if not ok: break
        if ok: V.append(r)
    vid = {r: i for i, r in enumerate(V)}

    # 2) cláusulas entre órbitas candidatas: recorrer vértices en el cuadrante
    allc = F + [c for r in V for c in cells(r)]
    owner = [-1] * len(F) + [vid[r] for r in V for _ in range(4)]
    A = np.array(allc); own = np.array(owner)
    clauses = set()
    for ia, a in enumerate(allc):
        if not (a[0] < h and a[1] < h): continue
        d = ((A - np.array(a)) ** 2).sum(1); d[ia] = -1
        order = np.argsort(d, kind='stable'); ds = d[order]
        cuts = np.flatnonzero(np.diff(ds)) + 1
        for grp in np.split(order, cuts):
            if len(grp) < 2 or d[grp[0]] <= 0: continue
            for i in range(len(grp)):
                for j in range(i + 1, len(grp)):
                    s = {own[ia], own[grp[i]], own[grp[j]]} - {-1}
                    if s: clauses.add(frozenset(s))
    unit = {next(iter(c)) for c in clauses if len(c) == 1}
    clauses = [c for c in clauses if len(c) > 1 and not (c & unit)]
    t1 = time.time()

    mdl = cp_model.CpModel()
    x = [mdl.NewBoolVar(f"o{i}") for i in range(len(V))]
    for u in unit: mdl.Add(x[u] == 0)
    for c in clauses: mdl.AddBoolOr([x[i].Not() for i in c])
    mdl.Maximize(sum(x))
    for r in V:                                   # pista: la solución actual
        mdl.AddHint(x[vid[r]], 1 if r in P else 0)
    sol = cp_model.CpSolver(); sol.parameters.max_time_in_seconds = secs; sol.parameters.num_workers = 8
    st = sol.Solve(mdl)
    best = int(sol.ObjectiveValue()); bound = int(sol.BestObjectiveBound() + 1e-9)
    estado = {cp_model.OPTIMAL: "ÓPTIMO DEMOSTRADO", cp_model.FEASIBLE: "factible (sin demostrar)"}.get(st, sol.StatusName(st))
    total = len(F) + 4 * best
    print(f"región {sys.argv[3]}: {len(R)} órbitas, {len(V)} candidatas, {len(clauses)} cláusulas "
          f"(armado {t1 - t0:.0f}s) | antes {old} órbitas, ahora {best} (cota {bound}) {estado} "
          f"| total {len(F) + 4 * old} -> {total}{'  *** MÁS PUNTOS ***' if best > old else ''}", flush=True)
    if out and st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        with open(out, 'w') as f:
            for p in F: f.write(f"{p[0]} {p[1]}\n")
            for r in V:
                if sol.Value(x[vid[r]]):
                    for c in cells(r): f.write(f"{c[0]} {c[1]}\n")


if __name__ == "__main__":
    main()
