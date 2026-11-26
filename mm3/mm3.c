/* mm3.c — n puntos en el espacio minimizando r^2 = (dist. máxima / dist. mínima)^2.
 *
 * Óptimo local: se normaliza para que la distancia máxima sea 1 y se busca
 * subir el umbral t de la distancia mínima al cuadrado. Para cada t de prueba se
 * minimiza con L-BFGS la penalización
 *     E = sum_{d^2 < t} (t - d^2)^2 + sum_{d^2 > 1} (d^2 - 1)^2
 * y se acepta el resultado sólo si su cociente REAL (calculado de las
 * coordenadas) es mejor; el paso en t se dobla o se parte a la mitad.
 * Búsqueda global: saltos entre cuencas (perturbaciones, reubicaciones y
 * arranques nuevos) con aceptación tipo Metropolis.
 *
 * Uso: ./mm3 n segundos semilla [inicial.txt|-] [récord_r2]
 * Escribe best_n<n>_seed<semilla>.txt cada vez que mejora.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAXN 128
#define MAXV (3 * MAXN)
#define M_LB 10

static int n, nv;
static unsigned long long rng;

static inline unsigned long long xr(void) {
    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng;
}
static inline double urand(void) { return (xr() >> 11) * (1.0 / 9007199254740992.0); }
static double nrand(void) {
    double u = urand() + 1e-300, v = urand();
    return sqrt(-2 * log(u)) * cos(2 * M_PI * v);
}
static double now(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* cociente real r^2 de una configuración; devuelve también min y max d^2 */
static double ratio(const double *X, double *dmn, double *dmx) {
    double mn = 1e300, mx = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            double a = X[3*i] - X[3*j], b = X[3*i+1] - X[3*j+1], c = X[3*i+2] - X[3*j+2];
            double d = a*a + b*b + c*c;
            if (d < mn) mn = d;
            if (d > mx) mx = d;
        }
    if (dmn) *dmn = mn;
    if (dmx) *dmx = mx;
    return mn > 0 ? mx / mn : 1e300;
}

/* centra y escala para que la distancia máxima sea 1 */
static void normalize(double *X) {
    double cx = 0, cy = 0, cz = 0, mx;
    for (int i = 0; i < n; i++) { cx += X[3*i]; cy += X[3*i+1]; cz += X[3*i+2]; }
    cx /= n; cy /= n; cz /= n;
    for (int i = 0; i < n; i++) { X[3*i] -= cx; X[3*i+1] -= cy; X[3*i+2] -= cz; }
    ratio(X, NULL, &mx);
    double s = 1.0 / sqrt(mx);
    for (int k = 0; k < nv; k++) X[k] *= s;
}

static double T_TGT;
static double energy(const double *X, double *G) {
    double E = 0, t = T_TGT;
    memset(G, 0, sizeof(double) * nv);
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            double a = X[3*i] - X[3*j], b = X[3*i+1] - X[3*j+1], c = X[3*i+2] - X[3*j+2];
            double d = a*a + b*b + c*c, f;
            if (d < t) { double v = t - d; E += v * v; f = -4 * v; }
            else if (d > 1) { double w = d - 1; E += w * w; f = 4 * w; }
            else continue;
            G[3*i] += f * a; G[3*i+1] += f * b; G[3*i+2] += f * c;
            G[3*j] -= f * a; G[3*j+1] -= f * b; G[3*j+2] -= f * c;
        }
    return E;
}

/* L-BFGS con búsqueda lineal de retroceso (Armijo) */
static double lbfgs(double *X, int maxit) {
    static double G[MAXV], Gn[MAXV], D[MAXV], Xn[MAXV];
    static double S[M_LB][MAXV], Y[M_LB][MAXV], rho[M_LB], al[M_LB];
    int k = 0, head = 0;
    double E = energy(X, G);
    for (int it = 0; it < maxit && E > 1e-32; it++) {
        /* dirección: -H g */
        memcpy(D, G, sizeof(double) * nv);
        int m = k < M_LB ? k : M_LB;
        for (int q = 0; q < m; q++) {
            int idx = (head - 1 - q + M_LB) % M_LB;
            double s = 0; for (int v = 0; v < nv; v++) s += S[idx][v] * D[v];
            al[idx] = rho[idx] * s;
            for (int v = 0; v < nv; v++) D[v] -= al[idx] * Y[idx][v];
        }
        double gam = 1;
        if (m > 0) {
            int idx = (head - 1 + M_LB) % M_LB;
            double sy = 0, yy = 0;
            for (int v = 0; v < nv; v++) { sy += S[idx][v] * Y[idx][v]; yy += Y[idx][v] * Y[idx][v]; }
            gam = sy / yy;
        } else {
            double gn = 0; for (int v = 0; v < nv; v++) gn += G[v] * G[v];
            gam = 1e-3 / (sqrt(gn) + 1e-300);
        }
        for (int v = 0; v < nv; v++) D[v] *= gam;
        for (int q = m - 1; q >= 0; q--) {
            int idx = (head - 1 - q + M_LB) % M_LB;
            double s = 0; for (int v = 0; v < nv; v++) s += Y[idx][v] * D[v];
            double be = rho[idx] * s;
            for (int v = 0; v < nv; v++) D[v] += S[idx][v] * (al[idx] - be);
        }
        double gd = 0; for (int v = 0; v < nv; v++) gd += G[v] * D[v];
        if (gd <= 0) { k = 0; head = 0; gd = 0; for (int v = 0; v < nv; v++) { D[v] = G[v] * 1e-3; gd += G[v] * D[v]; } }
        double step = 1, En = 0;
        int ok = 0;
        for (int ls = 0; ls < 40; ls++) {
            for (int v = 0; v < nv; v++) Xn[v] = X[v] - step * D[v];
            En = energy(Xn, Gn);
            if (En <= E - 1e-4 * step * gd) { ok = 1; break; }
            step *= 0.5;
        }
        if (!ok) break;
        double sy = 0;
        for (int v = 0; v < nv; v++) {
            S[head][v] = Xn[v] - X[v]; Y[head][v] = Gn[v] - G[v]; sy += S[head][v] * Y[head][v];
        }
        if (sy > 1e-300) { rho[head] = 1 / sy; head = (head + 1) % M_LB; k++; }
        memcpy(X, Xn, sizeof(double) * nv); memcpy(G, Gn, sizeof(double) * nv);
        if (E - En < 1e-18 * E && it > 50) { E = En; break; }
        E = En;
    }
    return E;
}

/* óptimo local: sube el umbral t de la distancia mínima; devuelve el r^2 real */
static double local_opt(double *X) {
    static double Yc[MAXV];
    normalize(X);
    double mn, r = ratio(X, &mn, NULL);
    double step = 0.05 * mn;
    for (int iter = 0; iter < 400 && step > 1e-15 * mn; iter++) {
        memcpy(Yc, X, sizeof(double) * nv);
        T_TGT = mn + step;
        lbfgs(Yc, 3000);
        normalize(Yc);
        double mn2, r2 = ratio(Yc, &mn2, NULL);
        if (r2 < r) { memcpy(X, Yc, sizeof(double) * nv); r = r2; mn = mn2; step *= 2; }
        else step *= 0.5;
    }
    return r;
}

static void random_config(double *X) {
    for (int i = 0; i < n; i++) {
        double x, y, z;
        do { x = 2 * urand() - 1; y = 2 * urand() - 1; z = 2 * urand() - 1; } while (x*x + y*y + z*z > 1);
        X[3*i] = 0.5 * x; X[3*i+1] = 0.5 * y; X[3*i+2] = 0.5 * z;
    }
}

static void save(const char *fn, const double *X, double r) {
    FILE *f = fopen(fn, "w");
    fprintf(f, "# n=%d r2=%.17g\n", n, r);
    for (int i = 0; i < n; i++) fprintf(f, "%.17g %.17g %.17g\n", X[3*i], X[3*i+1], X[3*i+2]);
    fclose(f);
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "uso: %s n segundos semilla [inicial|-] [récord_r2]\n", argv[0]); return 1; }
    n = atoi(argv[1]); nv = 3 * n;
    double secs = atof(argv[2]);
    unsigned long long seed = strtoull(argv[3], NULL, 10);
    rng = seed * 0x9E3779B97F4A7C15ULL + 0x1234567ULL; if (!rng) rng = 1;
    for (int i = 0; i < 20; i++) xr();
    const char *init = argc > 4 && strcmp(argv[4], "-") ? argv[4] : NULL;
    double rec = argc > 5 ? atof(argv[5]) : 0;
    static double cur[MAXV], best[MAXV], trial[MAXV];
    char fn[256]; snprintf(fn, sizeof fn, "best_n%d_seed%llu.txt", n, seed);

    if (init) {
        FILE *f = fopen(init, "r"); char line[512]; int i = 0;
        while (i < n && fgets(line, sizeof line, f)) {
            if (line[0] == '#') continue;
            if (sscanf(line, "%lf %lf %lf", &cur[3*i], &cur[3*i+1], &cur[3*i+2]) == 3) i++;
        }
        fclose(f);
        if (i != n) { fprintf(stderr, "inicial: leí %d puntos, esperaba %d\n", i, n); return 1; }
    } else random_config(cur);
    double rc = local_opt(cur), rb = rc;
    memcpy(best, cur, sizeof(double) * nv);
    save(fn, best, rb);
    printf("inicio r2=%.12f\n", rb); fflush(stdout);

    double t0 = now(), Temp = 1e-3 * rb;
    long long hops = 0, since = 0;
    while (now() - t0 < secs) {
        memcpy(trial, cur, sizeof(double) * nv);
        double u = urand(), mn;
        ratio(trial, &mn, NULL);
        double sc = sqrt(mn);
        if (u < 0.45) {                       /* sacudir todo */
            double a = sc * (0.02 + 0.25 * urand());
            for (int v = 0; v < nv; v++) trial[v] += a * nrand();
        } else if (u < 0.8) {                 /* reubicar 1–3 puntos */
            int c = 1 + (int)(3 * urand());
            for (int q = 0; q < c; q++) {
                int i = (int)(n * urand());
                double x, y, z;
                do { x = 2 * urand() - 1; y = 2 * urand() - 1; z = 2 * urand() - 1; } while (x*x + y*y + z*z > 1);
                trial[3*i] = 0.5 * x; trial[3*i+1] = 0.5 * y; trial[3*i+2] = 0.5 * z;
            }
        } else if (u < 0.9) {                 /* mover 1 punto un poco */
            int i = (int)(n * urand());
            for (int k = 0; k < 3; k++) trial[3*i+k] += sc * 0.5 * nrand();
        } else random_config(trial);          /* arranque nuevo */
        double rt = local_opt(trial);
        hops++; since++;
        if (rt < rc || urand() < exp(-(rt - rc) / Temp)) { memcpy(cur, trial, sizeof(double) * nv); rc = rt; }
        if (rt < rb - 1e-13) {
            rb = rt; memcpy(best, trial, sizeof(double) * nv); save(fn, best, rb); since = 0;
            printf("%.0fs salto %lld: r2=%.13f%s\n", now() - t0, hops, rb,
                   rec > 0 && rb < rec ? "  < récord" : ""); fflush(stdout);
        }
        if (since > 0 && since % 2000 == 0) { memcpy(cur, best, sizeof(double) * nv); rc = rb; }
    }
    printf("FIN n=%d semilla=%llu mejor r2=%.15f (saltos %lld)\n", n, seed, rb, hops);
    return 0;
}
