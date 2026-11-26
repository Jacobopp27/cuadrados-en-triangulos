// isos.c — búsqueda local (tabú) del mayor subconjunto de la cuadrícula n×n
// sin triángulos isósceles (problema 6.59 de AlphaEvolve).
//
// Violaciones V(S) = número de pares (ápice a ∈ S, par {b,c} ⊂ S) con
// d(a,b) = d(a,c). S es válido si V = 0. Para cada celda x se mantienen:
//   A[x] = pares {b,c} ⊂ S equidistantes de x   (x sería el ápice)
//   B[x] = pares (a,c) de S con d(a,c) = d(a,x)  (x sería parte de la base)
// Se actualizan recorriendo mediatrices y círculos, sin tablas grandes.
// Añadir x cuesta A[x] + B[x] violaciones nuevas. Con |S| = k fijo se minimiza
// V intercambiando puntos; al llegar a V = 0 se guarda y se pasa a k + 1.
//
// Uso: ./isos n segundos semilla [inicial.txt] [simetria]
//   simetria: 0 ninguna, 1 central (giro 180°), 2 espejo x, 3 dos espejos,
//   (los giros de 90° no sirven: cada órbita sería un cuadrado, que es isósceles)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>

static int n, N, R, SYM;
static int *A, *B, *inS, *pos, *S, k;
static long long V;
static int *offs_start, *offs;  // offsets por radio: offs[2*i], offs[2*i+1]
static int *tabu_add, *tabu_rem;
static int *orb, *osz;          // órbita de cada celda bajo el grupo de simetría (hasta 8)

static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline int rndi(int m) { return (int)(rnd() % (uint64_t)m); }
static double now(void) { struct timeval tv; gettimeofday(&tv, NULL); return tv.tv_sec + 1e-6 * tv.tv_usec; }

static inline int d2(int a, int b) { int dx = a / n - b / n, dy = a % n - b % n; return dx * dx + dy * dy; }

// Aplica (signo +1/-1) los términos de los pares {p, b}, b ∈ S, b ≠ p:
//  (1) A[x] para x en la mediatriz de p y b
//  (2) B[x] para x en el círculo de centro p que pasa por b (x ≠ b)
//  (3) B[x] para x en el círculo de centro b que pasa por p (x ≠ p)
static inline void circle(int cx, int cy, int r, int skip, int sg) {
    for (int q = offs_start[r]; q < offs_start[r + 1]; q++) {
        int x = cx + offs[2 * q], y = cy + offs[2 * q + 1];
        if ((unsigned)x >= (unsigned)n || (unsigned)y >= (unsigned)n) continue;
        int c = x * n + y; if (c != skip) B[c] += sg;
    }
}
static void apply_pairs(int p, int sg) {
    int px = p / n, py = p % n;
    for (int t = 0; t < k; t++) {
        int b = S[t]; if (b == p) continue;
        int bx = b / n, by = b % n, u = bx - px, v = by - py, r = u * u + v * v;
        circle(px, py, r, b, sg);
        circle(bx, by, r, p, sg);
        int c = bx * bx + by * by - px * px - py * py;
        if (c & 1) continue;
        int c2 = c / 2;                       // u X + v Y = c2
        if (abs(v) >= abs(u)) {
            for (int X = 0; X < n; X++) { int w = c2 - u * X; if (w % v == 0) { int Y = w / v; if ((unsigned)Y < (unsigned)n) A[X * n + Y] += sg; } }
        } else {
            for (int Y = 0; Y < n; Y++) { int w = c2 - v * Y; if (w % u == 0) { int X = w / u; if ((unsigned)X < (unsigned)n) A[X * n + Y] += sg; } }
        }
    }
}

static void add_pt(int p) {
    V += A[p] + B[p];
    apply_pairs(p, +1);
    inS[p] = 1; pos[p] = k; S[k++] = p;
}

static void rem_pt(int p) {
    int i = pos[p]; S[i] = S[--k]; pos[S[i]] = i; inS[p] = 0;
    apply_pairs(p, -1);
    V -= A[p] + B[p];
}
static long long brute_V(void) {
    long long v = 0;
    for (int t = 0; t < k; t++) { int a = S[t];
        for (int u = 0; u < k; u++) for (int w = u + 1; w < k; w++) if (S[u] != a && S[w] != a && d2(a, S[u]) == d2(a, S[w])) v++; }
    return v;
}

// costo de quitar u ∈ S
static inline long long rem_score(int u) { return (long long)A[u] + B[u]; }

static void save(const char *fn) {
    FILE *f = fopen(fn, "w"); if (!f) return;
    fprintf(f, "# n=%d puntos=%d\n", n, k);
    for (int t = 0; t < k; t++) fprintf(f, "%d %d\n", S[t] / n, S[t] % n);
    fclose(f);
}

// imagen de la celda (x,y) bajo el elemento g del grupo de simetría del cuadrado
static int img(int g, int x, int y) {
    int m = n - 1, X, Y;
    switch (g) {
        case 0: X = x; Y = y; break;
        case 1: X = m - x; Y = m - y; break;      // giro 180°
        case 2: X = m - x; Y = y; break;          // espejo x
        case 3: X = x; Y = m - y; break;          // espejo y
        case 4: X = m - y; Y = x; break;          // giro 90°
        case 5: X = y; Y = m - x; break;          // giro 270°
        case 6: X = y; Y = x; break;              // diagonal
        default: X = m - y; Y = m - x; break;     // antidiagonal
    }
    return X * n + Y;
}

// costo aproximado de añadir la órbita de c (suma de costos individuales)
static long long orbit_cost(int c, long long it, int check_tabu) {
    long long cc = 0;
    for (int j = 0; j < osz[c]; j++) {
        int e = orb[8 * c + j];
        if (inS[e]) return -1;
        if (check_tabu && tabu_add[e] > it) return -1;
        cc += A[e] + B[e];
    }
    return cc;
}
static void add_orbit(int c) { for (int j = 0; j < osz[c]; j++) { int e = orb[8 * c + j]; if (!inS[e]) add_pt(e); } }
static void rem_orbit(int c) { for (int j = 0; j < osz[c]; j++) { int e = orb[8 * c + j]; if (inS[e]) rem_pt(e); } }

static int best_add(long long it, int check_tabu) {
    int bx = -1, ties = 0; long long bc = 1LL << 60;
    for (int c = 0; c < N; c++) {
        if (inS[c] || orb[8 * c] != c) continue;   // solo el representante de cada órbita
        long long cc = orbit_cost(c, it, check_tabu);
        if (cc < 0) continue;
        cc = cc * 8 / osz[c];                         // por punto, para comparar órbitas de distinto tamaño
        if (cc < bc) { bc = cc; bx = c; ties = 1; } else if (cc == bc && rndi(++ties) == 0) bx = c;
    }
    return bx;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "uso: %s n segundos semilla [inicial.txt|-] [simetria]\n", argv[0]); return 1; }
    n = atoi(argv[1]); double budget = atof(argv[2]); unsigned long long seed = strtoull(argv[3], NULL, 10);
    const char *init = argc > 4 && strcmp(argv[4], "-") ? argv[4] : NULL;
    SYM = argc > 5 ? atoi(argv[5]) : 0;
    rs = 0x9E3779B97F4A7C15ULL ^ (seed * 0xBF58476D1CE4E5B9ULL); for (int i = 0; i < 20; i++) rnd();
    N = n * n; R = 2 * (n - 1) * (n - 1) + 1;
    A = calloc(N, sizeof *A); B = calloc(N, sizeof *B); inS = calloc(N, sizeof *inS);
    pos = calloc(N, sizeof *pos); S = calloc(N, sizeof *S);
    tabu_add = calloc(N, sizeof *tabu_add); tabu_rem = calloc(N, sizeof *tabu_rem);
    orb = malloc(8 * N * sizeof *orb); osz = malloc(N * sizeof *osz);
    if (!A || !B || !orb) { fprintf(stderr, "sin memoria\n"); return 1; }
    static const int groups[6][8] = {{0}, {0, 1}, {0, 2}, {0, 1, 2, 3}, {0, 1, 4, 5}, {0, 1, 2, 3, 4, 5, 6, 7}};
    static const int gsize[6] = {1, 2, 2, 4, 4, 8};
    if (SYM < 0 || SYM > 5) SYM = 0;
    for (int c = 0; c < N; c++) {
        int x = c / n, y = c % n, m = 0, tmp[8];
        for (int g = 0; g < gsize[SYM]; g++) {
            int e = img(groups[SYM][g], x, y), dup = 0;
            for (int j = 0; j < m; j++) if (tmp[j] == e) dup = 1;
            if (!dup) tmp[m++] = e;
        }
        // representante = menor índice, en la primera posición
        for (int j = 1; j < m; j++) if (tmp[j] < tmp[0]) { int s = tmp[0]; tmp[0] = tmp[j]; tmp[j] = s; }
        for (int j = 0; j < m; j++) orb[8 * c + j] = tmp[j];
        osz[c] = m;
    }
    // offsets por radio
    offs_start = calloc(R + 2, sizeof *offs_start);
    for (int dx = -(n - 1); dx <= n - 1; dx++) for (int dy = -(n - 1); dy <= n - 1; dy++) offs_start[dx * dx + dy * dy + 1]++;
    for (int r = 0; r <= R; r++) offs_start[r + 1] += offs_start[r];
    offs = malloc(2 * sizeof *offs * (size_t)offs_start[R + 1]);
    int *fillp = calloc(R + 1, sizeof *fillp);
    for (int dx = -(n - 1); dx <= n - 1; dx++) for (int dy = -(n - 1); dy <= n - 1; dy++) {
        int r = dx * dx + dy * dy, q = offs_start[r] + fillp[r]++;
        offs[2 * q] = dx; offs[2 * q + 1] = dy;
    }

    char out[256]; snprintf(out, sizeof out, "best_n%d_sym%d_seed%llu.txt", n, SYM, seed);
    if (init) {
        FILE *f = fopen(init, "r"); char line[256]; int x, y;
        while (f && fgets(line, sizeof line, f)) if (line[0] != '#' && sscanf(line, "%d %d", &x, &y) == 2 && !inS[x * n + y]) add_pt(x * n + y);
        if (f) fclose(f);
        printf("inicial: %d puntos, V = %lld\n", k, V);
    }
    int best = V == 0 ? k : 0;
    double t0 = now(); long long it = 0, last_imp = 0; long long bestV = V;
    int tt_add = 10, tt_rem = 3;
    while (now() - t0 < budget) {
        for (int rep = 0; rep < 256; rep++) {
            it++;
            if (V == 0) {
                if (k > best) { best = k; save(out); printf("[t=%.0fs it=%lld] %d puntos válidos\n", now() - t0, it, k); fflush(stdout); }
                int bx = best_add(it, 0);
                if (bx < 0) break;
                add_orbit(bx);
                bestV = V; last_imp = it;
                continue;
            }
            // quitar: la órbita del punto de mayor costo no tabú (a veces uno al azar en conflicto)
            int u = -1;
            if (rndi(100) < 70) {
                long long bs = -1; int ties = 0;
                for (int t = 0; t < k; t++) {
                    int c = S[t]; if (tabu_rem[c] > it) continue;
                    long long s = rem_score(c);
                    if (s > bs) { bs = s; u = c; ties = 1; } else if (s == bs && rndi(++ties) == 0) u = c;
                }
            }
            if (u < 0) { for (int tries = 0; tries < 50; tries++) { int c = S[rndi(k)]; if (rem_score(c) > 0) { u = c; break; } } }
            if (u < 0) u = S[rndi(k)];
            int tu = it + tt_add + rndi(10);
            for (int j = 0; j < osz[u]; j++) tabu_add[orb[8 * u + j]] = tu;
            rem_orbit(u);
            int bx = best_add(it, 1);
            if (bx < 0) continue;
            add_orbit(bx);
            int tr = it + tt_rem + rndi(3);
            for (int j = 0; j < osz[bx]; j++) tabu_rem[orb[8 * bx + j]] = tr;
            if (V < bestV) { bestV = V; last_imp = it; }
#ifdef CHECK
            { long long bv = brute_V(); if (bv != V) { printf("ERROR it=%lld V=%lld brute=%lld k=%d\n", it, V, bv, k); return 1; } }
#endif
            // estancamiento: perturbación (quitar varias órbitas al azar y rellenar con las más baratas)
            if (it - last_imp > 20000 + 200LL * k) {
                int kk = k, m = 2 + rndi(4);
                for (int q = 0; q < m && k > 1; q++) { int c = S[rndi(k)]; for (int j = 0; j < osz[c]; j++) tabu_add[orb[8 * c + j]] = it + tt_add; rem_orbit(c); }
                while (k < kk) { int bb = best_add(it, 1); if (bb < 0) break; add_orbit(bb); }
                bestV = V; last_imp = it;
            }
        }
    }
    printf("FIN n=%d semilla=%llu simetria=%d mejor = %d puntos (iteraciones %lld, %.0f it/s)\n", n, seed, SYM, best, it, it / (now() - t0));
    return 0;
}
