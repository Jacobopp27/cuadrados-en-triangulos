"""exacto.py — valor EXACTO (con demostración) del mayor conjunto sin isósceles en la
cuadrícula n×n, sin suponer simetría. Usa CP-SAT (OR-Tools).

Restricción: para cada vértice a y cada par {b, c} con d(a,b) = d(a,c), no pueden
estar los tres (ver abajo el caso degenerado). Se rompe parte de la
simetría del cuadrado (8 rotaciones/reflexiones) con restricciones lexicográficas
sobre los conteos de puntos por fila/columna extremas.

Versión OEIS A271907 (por defecto): NO se cuenta el caso degenerado (a punto medio de b y c),
porque tres puntos en línea no forman triángulo. Con --degenerado se cuenta (versión AlphaEvolve 6.59).

Uso: python exacto.py n segundos [hilos] [solucion_salida.txt] [--degenerado]
"""
import sys, time
from collections import defaultdict
from ortools.sat.python import cp_model


def main():
    deg = '--degenerado' in sys.argv; sys.argv = [a for a in sys.argv if a != '--degenerado']
    n = int(sys.argv[1]); secs = float(sys.argv[2])
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    out = sys.argv[4] if len(sys.argv) > 4 else None
    cells = [(x, y) for x in range(n) for y in range(n)]
    idx = {c: i for i, c in enumerate(cells)}
    t0 = time.time()
    tri = set()
    for a in cells:
        groups = defaultdict(list)
        for b in cells:
            if b != a:
                groups[(a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2].append(idx[b])
        ia = idx[a]
        for g in groups.values():
            for i in range(len(g)):
                for j in range(i + 1, len(g)):
                    b, c = cells[g[i]], cells[g[j]]
                    if not deg and b[0] + c[0] == 2 * a[0] and b[1] + c[1] == 2 * a[1]:
                        continue   # a es punto medio de b y c: colineales, no es triángulo
                    tri.add(tuple(sorted((ia, g[i], g[j]))))
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"c{i}") for i in range(len(cells))]
    for t in tri:
        m.AddBoolOr([x[i].Not() for i in t])
    tot = sum(x)
    # ruptura de simetría: la fila 0 tiene al menos tantos puntos como las otras 3 líneas
    # del borde, y la columna 0 al menos tantos como la columna n-1
    row0 = sum(x[idx[(i, 0)]] for i in range(n)); rowN = sum(x[idx[(i, n - 1)]] for i in range(n))
    col0 = sum(x[idx[(0, i)]] for i in range(n)); colN = sum(x[idx[(n - 1, i)]] for i in range(n))
    m.Add(row0 >= rowN); m.Add(row0 >= col0); m.Add(row0 >= colN); m.Add(col0 >= colN)
    m.Maximize(tot)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = secs; s.parameters.num_workers = workers
    st = s.Solve(m)
    estado = {cp_model.OPTIMAL: "ÓPTIMO DEMOSTRADO", cp_model.FEASIBLE: "factible (sin demostrar)"}.get(st, s.StatusName(st))
    print(f"n={n} ({'con' if deg else 'sin'} degenerado): {len(tri)} tríos prohibidos | mejor {int(s.ObjectiveValue())}, cota {int(s.BestObjectiveBound() + 1e-9)} "
          f"| {estado} | {time.time() - t0:.0f}s", flush=True)
    if out and st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        with open(out, "w") as f:
            for i, c in enumerate(cells):
                if s.Value(x[i]): f.write(f"{c[0]} {c[1]}\n")


if __name__ == "__main__":
    main()
