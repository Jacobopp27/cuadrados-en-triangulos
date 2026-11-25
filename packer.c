// packer.c — busca empaquetamientos de n cuadrados unitarios dentro del
// triángulo equilátero más pequeño posible (problema "Squares in Triangles"
// del Packing Center de Erich Friedman).
//
// Método: energía de penalización (solapamiento entre cuadrados + salida del
// triángulo), minimizada con FIRE; se comprime el triángulo mientras haya una
// configuración factible y se reinicia desde configuraciones aleatorias.
//
// Uso: ./packer n segundos semilla s_inicial [archivo_inicial]
// Escribe la mejor configuración en best_n<n>_seed<semilla>.txt

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 128
#ifndef MAXIT
#define MAXIT 10000
#endif
static const double R3 = 1.7320508075688772935;

static int N;
static double X[MAXN], Y[MAXN], T[MAXN];
static double CS[MAXN], SN[MAXN];
static double S; // lado del triángulo

static unsigned long long rng_state;
static double urand(void) {
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17;
    return (rng_state >> 11) * (1.0 / 9007199254740992.0);
}

static double now(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// Profundidad de penetración (SAT) entre dos cuadrados unitarios; 0 si no se tocan.
static double pair_depth(double xi, double yi, double ci, double si,
                         double xj, double yj, double cj, double sj) {
    double dx = xj - xi, dy = yj - yi;
    if (dx * dx + dy * dy >= 2.0) return 0.0;
    double ax[4][2] = {{ci, si}, {-si, ci}, {cj, sj}, {-sj, cj}};
    double best = 1e300;
    for (int k = 0; k < 4; k++) {
        double a0 = ax[k][0], a1 = ax[k][1];
        double ri = 0.5 * (fabs(a0 * ci + a1 * si) + fabs(-a0 * si + a1 * ci));
        double rj = 0.5 * (fabs(a0 * cj + a1 * sj) + fabs(-a0 * sj + a1 * cj));
        double o = ri + rj - fabs(a0 * dx + a1 * dy);
        if (o <= 0) return 0.0;
        if (o < best) best = o;
    }
    return best;
}

// Violación de contención de un cuadrado (suma de cuadrados de lo que sale por cada lado).
static double cont_energy(double x, double y, double c, double s) {
    // normales hacia adentro y offsets: n·p >= off
    double nx[3] = {0.0, R3 / 2, -R3 / 2};
    double ny[3] = {1.0, -0.5, -0.5};
    double off[3] = {0.0, 0.0, -R3 / 2 * S};
    double e = 0;
    for (int k = 0; k < 3; k++) {
        double r = 0.5 * (fabs(nx[k] * c + ny[k] * s) + fabs(-nx[k] * s + ny[k] * c));
        double v = off[k] - (nx[k] * x + ny[k] * y - r);
        if (v > 0) e += v * v;
    }
    return e;
}

static double local_energy(int i, double x, double y, double t) {
    double c = cos(t), s = sin(t);
    double e = cont_energy(x, y, c, s);
    for (int j = 0; j < N; j++) {
        if (j == i) continue;
        double d = pair_depth(x, y, c, s, X[j], Y[j], CS[j], SN[j]);
        e += d * d;
    }
    return e;
}

static void update_trig(void) { for (int i = 0; i < N; i++) { CS[i] = cos(T[i]); SN[i] = sin(T[i]); } }

static double total_energy(double *maxdepth) {
    double e = 0, md = 0;
    for (int i = 0; i < N; i++) {
        double ce = cont_energy(X[i], Y[i], CS[i], SN[i]);
        e += ce; if (sqrt(ce) > md) md = sqrt(ce);
        for (int j = i + 1; j < N; j++) {
            double d = pair_depth(X[i], Y[i], CS[i], SN[i], X[j], Y[j], CS[j], SN[j]);
            e += d * d; if (d > md) md = d;
        }
    }
    if (maxdepth) *maxdepth = md;
    return e;
}

static void gradient_fd(double *g) {
    const double h = 1e-7;
    for (int i = 0; i < N; i++) {
        double x = X[i], y = Y[i], t = T[i];
        g[3 * i + 0] = (local_energy(i, x + h, y, t) - local_energy(i, x - h, y, t)) / (2 * h);
        g[3 * i + 1] = (local_energy(i, x, y + h, t) - local_energy(i, x, y - h, t)) / (2 * h);
        g[3 * i + 2] = (local_energy(i, x, y, t + h) - local_energy(i, x, y, t - h)) / (2 * h);
    }
}

static inline double sgn(double v) { return (v > 0) - (v < 0); }

// Gradiente analítico de la energía total (requiere CS/SN actualizados).
static void gradient(double *g) {
    memset(g, 0, sizeof(double) * 3 * N);
    const double nx[3] = {0.0, R3 / 2, -R3 / 2};
    const double ny[3] = {1.0, -0.5, -0.5};
    const double off[3] = {0.0, 0.0, -R3 / 2 * S};
    for (int i = 0; i < N; i++) {
        double c = CS[i], s = SN[i];
        for (int k = 0; k < 3; k++) {
            double nu = nx[k] * c + ny[k] * s, nv = -nx[k] * s + ny[k] * c;
            double r = 0.5 * (fabs(nu) + fabs(nv));
            double v = off[k] - (nx[k] * X[i] + ny[k] * Y[i] - r);
            if (v > 0) {
                g[3 * i] += 2 * v * (-nx[k]);
                g[3 * i + 1] += 2 * v * (-ny[k]);
                g[3 * i + 2] += 2 * v * 0.5 * (sgn(nu) * nv - sgn(nv) * nu);
            }
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dx = X[j] - X[i], dy = Y[j] - Y[i];
            if (dx * dx + dy * dy >= 2.0) continue;
            // ejes: u_i, v_i, u_j, v_j
            double ax[4] = {CS[i], -SN[i], CS[j], -SN[j]};
            double ay[4] = {SN[i], CS[i], SN[j], CS[j]};
            double best = 1e300; int kb = -1; double pb = 0;
            int sep = 0;
            for (int k = 0; k < 4; k++) {
                double a0 = ax[k], a1 = ay[k];
                double ri = 0.5 * (fabs(a0 * CS[i] + a1 * SN[i]) + fabs(-a0 * SN[i] + a1 * CS[i]));
                double rj = 0.5 * (fabs(a0 * CS[j] + a1 * SN[j]) + fabs(-a0 * SN[j] + a1 * CS[j]));
                double p = a0 * dx + a1 * dy;
                double o = ri + rj - fabs(p);
                if (o <= 0) { sep = 1; break; }
                if (o < best) { best = o; kb = k; pb = p; }
            }
            if (sep) continue;
            double d = best, a0 = ax[kb], a1 = ay[kb];
            int owner = kb < 2 ? i : j, m = kb < 2 ? j : i;
            double cphi = a0 * CS[m] + a1 * SN[m];
            double sphi = -a0 * SN[m] + a1 * CS[m];
            double drdphi = 0.5 * (-sgn(cphi) * sphi + sgn(sphi) * cphi);
            double sp = sgn(pb);
            double adx = -a1 * dx + a0 * dy; // a'·Δ
            double f = 2 * d;
            g[3 * i] += f * sp * a0;  g[3 * i + 1] += f * sp * a1;
            g[3 * j] -= f * sp * a0;  g[3 * j + 1] -= f * sp * a1;
            g[3 * owner + 2] += f * (drdphi - sp * adx);
            g[3 * m + 2] += f * (-drdphi);
        }
    }
}

// Minimización FIRE. Devuelve la energía final.
static double minimize(int maxit, double tol) {
    static double v[3 * MAXN], g[3 * MAXN];
    int D = 3 * N;
    memset(v, 0, sizeof(double) * D);
    double dt = 0.02, dtmax = 0.1, alpha = 0.1;
    int npos = 0;
    double E = 0;
    for (int it = 0; it < maxit; it++) {
        update_trig();
        if ((it & 31) == 0) { E = total_energy(NULL); if (E < tol) return E; }
        gradient(g);
        double P = 0, vn = 0, fn = 0;
        for (int k = 0; k < D; k++) { P += -g[k] * v[k]; vn += v[k] * v[k]; fn += g[k] * g[k]; }
        vn = sqrt(vn); fn = sqrt(fn);
        if (fn < 1e-30) break;
        for (int k = 0; k < D; k++) v[k] = (1 - alpha) * v[k] + alpha * (-g[k]) / fn * vn;
        if (P > 0) {
            if (++npos > 5) { dt = fmin(dt * 1.1, dtmax); alpha *= 0.99; }
        } else {
            npos = 0; dt *= 0.5; alpha = 0.1;
            memset(v, 0, sizeof(double) * D);
        }
        for (int k = 0; k < D; k++) v[k] += -g[k] * dt;
        for (int i = 0; i < N; i++) { X[i] += v[3 * i] * dt; Y[i] += v[3 * i + 1] * dt; T[i] += v[3 * i + 2] * dt; }
    }
    update_trig();
    return total_energy(NULL);
}

// Lado mínimo del triángulo (con la misma orientación) que contiene todos los
// vértices; además traslada la configuración a la posición canónica.
static double fit_and_canon(void) {
    double minY = 1e300, minL = 1e300, maxR = -1e300;
    update_trig();
    for (int i = 0; i < N; i++) {
        double c = CS[i], s = SN[i];
        for (int a = -1; a <= 1; a += 2) for (int b = -1; b <= 1; b += 2) {
            double vx = X[i] + 0.5 * a * c - 0.5 * b * s;
            double vy = Y[i] + 0.5 * a * s + 0.5 * b * c;
            if (vy < minY) minY = vy;
            if (vx - vy / R3 < minL) minL = vx - vy / R3;
            if (vx + vy / R3 > maxR) maxR = vx + vy / R3;
        }
    }
    // traslación: y -> y - minY ; luego x tal que el lado izquierdo toque
    for (int i = 0; i < N; i++) Y[i] -= minY;
    double shiftx = -(minL + minY / R3);
    for (int i = 0; i < N; i++) X[i] += shiftx;
    return maxR - minL - 2 * minY / R3;
}

static void random_init(double s) {
    S = s;
    for (int i = 0; i < N; i++) {
        double u, w;
        do { u = urand(); w = urand(); } while (u + w > 1);
        // punto uniforme en el triángulo reducido
        double m = 0.7, sh = s - 2 * m * R3; if (sh < 0.1) sh = 0.1;
        double ox = m * R3, oy = m;
        X[i] = ox + sh * (u + 0.5 * w);
        Y[i] = oy + sh * (R3 / 2) * w;
        // la mayoría de cuadrados alineados con algún lado del triángulo
        // (0°, 30°, 60° módulo 90°), como en los récords conocidos
        double r = urand();
        if (r < 0.85) {
            // ángulo del lado más cercano: abajo 0°, izquierdo 60°, derecho 120° (≡30°)
            double db = Y[i];
            double dl = (R3 / 2) * X[i] - 0.5 * Y[i];
            double dr = -(R3 / 2) * (X[i] - s) - 0.5 * Y[i];
            double ang = (db <= dl && db <= dr) ? 0.0 : (dl <= dr ? M_PI / 3 : M_PI / 6);
            if (urand() < 0.15) ang = (int)(urand() * 3) * (M_PI / 6);
            T[i] = ang + 0.02 * (urand() - 0.5);
        } else T[i] = urand() * M_PI / 2;
    }
    update_trig();
}

static void scale_to(double snew, double sold) {
    double f = snew / sold;
    for (int i = 0; i < N; i++) { X[i] *= f; Y[i] *= f; }
    S = snew;
}

static void save(const char *fn, double s, double md) {
    FILE *f = fopen(fn, "w");
    if (!f) return;
    fprintf(f, "%d %.15f %.3e\n", N, s, md);
    for (int i = 0; i < N; i++) fprintf(f, "%.15f %.15f %.15f\n", X[i], Y[i], T[i]);
    fclose(f);
}

static double best = 1e300, t0, budget;
static unsigned long long seed;
static char out[256];
static double gx[MAXN], gy[MAXN], gt[MAXN];
static const double TOL = 1e-22;

// Comprime desde la configuración actual empezando con lado s_start.
// Devuelve el menor lado factible (configuración en gx/gy/gt) o -1.
static double compress(double s_start, int allow_grow) {
    S = s_start;
    double good_s = -1, step = 0.003;
    int fails = 0, grow = 0;
    while (step > 1e-8 && now() - t0 < budget) {
        double E = minimize(MAXIT, TOL);
        if (E < TOL) {
            double sn = fit_and_canon();
            for (int i = 0; i < N; i++) { gx[i] = X[i]; gy[i] = Y[i]; gt[i] = T[i]; }
            good_s = sn;
            if (sn < best - 1e-12) {
                double md; S = sn; update_trig(); total_energy(&md);
                best = sn; save(out, sn, md);
                printf("[n=%d seed=%llu t=%.0fs] mejor s = %.9f\n", N, seed, now() - t0, sn);
                fflush(stdout);
            }
            step = fmin(step * 1.5, 0.02);
            S = sn; scale_to(sn * (1 - step), sn);
            fails = 0;
        } else if (good_s < 0) {
            if (!allow_grow && ++grow > 2) return -1;
            scale_to(S * (allow_grow ? 1.05 : 1.002), S);
        } else if (++fails < 3) {
            for (int i = 0; i < N; i++) { X[i] += 0.03 * (urand() - 0.5); Y[i] += 0.03 * (urand() - 0.5); T[i] += 0.08 * (urand() - 0.5); }
        } else {
            for (int i = 0; i < N; i++) { X[i] = gx[i]; Y[i] = gy[i]; T[i] = gt[i]; }
            step *= 0.5; fails = 0;
            S = good_s; scale_to(good_s * (1 - step), good_s);
        }
    }
    return good_s;
}

static void random_point(double s, double *x, double *y) {
    double u, w;
    do { u = urand(); w = urand(); } while (u + w > 1);
    double m = 0.6, sh = s - 2 * m * R3; if (sh < 0.1) sh = 0.1;
    *x = m * R3 + sh * (u + 0.5 * w);
    *y = m + sh * (R3 / 2) * w;
}

static void perturb(double s) {
    double r = urand();
    if (r < 0.5) {
        // reubicar 1 a 3 cuadrados en lugares aleatorios
        int k = 1 + (int)(urand() * 3);
        for (int q = 0; q < k; q++) {
            int i = (int)(urand() * N);
            random_point(s, &X[i], &Y[i]);
            T[i] = (urand() < 0.7) ? (int)(urand() * 3) * (M_PI / 6) : urand() * M_PI / 2;
        }
    } else if (r < 0.8) {
        double a = 0.05 + 0.2 * urand();
        for (int i = 0; i < N; i++) { X[i] += a * (urand() - 0.5); Y[i] += a * (urand() - 0.5); T[i] += a * (urand() - 0.5); }
    } else {
        int k = 1 + (int)(urand() * 4);
        for (int q = 0; q < k; q++) { int i = (int)(urand() * N); T[i] += (urand() - 0.5) * 0.8; }
    }
    update_trig();
}

// Busca una configuración factible con el lado FIJO s (basin hopping sobre la energía).
// Si la encuentra, la comprime. Devuelve 1 si encontró algo factible.
static int load_init(const char *fn, int drop);
static const char *fixed_init = NULL;
static int fixed_search(double s, int iters) {
    static double bx[MAXN], by[MAXN], bt[MAXN];
    static int attempt = 0;
    if (fixed_init) {
        // cuadrado quitado distinto en cada intento
        if (!load_init(fixed_init, (int)((seed + attempt) % (unsigned long long)(N + 1)))) exit(1);
        double sf = S; scale_to(s, sf);
        for (int i = 0; i < N; i++) { X[i] += 0.01 * (urand() - 0.5); Y[i] += 0.01 * (urand() - 0.5); }
        update_trig();
    } else random_init(s);
    attempt++;
    double Ecur = minimize(MAXIT, TOL);
    for (int i = 0; i < N; i++) { bx[i] = X[i]; by[i] = Y[i]; bt[i] = T[i]; }
    for (int it = 0; it < iters && now() - t0 < budget; it++) {
        if (Ecur < TOL) {
            printf("[n=%d seed=%llu t=%.0fs] ¡factible con s fijo = %.9f! comprimiendo...\n", N, seed, now() - t0, s);
            fflush(stdout);
            update_trig();
            compress(s, 0);
            return 1;
        }
        for (int i = 0; i < N; i++) { X[i] = bx[i]; Y[i] = by[i]; T[i] = bt[i]; }
        update_trig();
        double r = urand();
        if (r < 0.6) {
            // reubicar el cuadrado con más energía (y a veces otro al azar)
            int worst = 0; double we = -1;
            for (int i = 0; i < N; i++) { double e = local_energy(i, X[i], Y[i], T[i]); if (e > we) { we = e; worst = i; } }
            random_point(s, &X[worst], &Y[worst]);
            T[worst] = (urand() < 0.6) ? (int)(urand() * 3) * (M_PI / 6) : urand() * M_PI / 2;
            if (urand() < 0.3) { int i = (int)(urand() * N); random_point(s, &X[i], &Y[i]); T[i] = urand() * M_PI / 2; }
        } else {
            double a = 0.03 + 0.25 * urand();
            for (int i = 0; i < N; i++) { X[i] += a * (urand() - 0.5); Y[i] += a * (urand() - 0.5); T[i] += a * (urand() - 0.5); }
        }
        S = s;
        double E = minimize(MAXIT, TOL);
        if (E < Ecur) {
            Ecur = E;
            for (int i = 0; i < N; i++) { bx[i] = X[i]; by[i] = Y[i]; bt[i] = T[i]; }
        }
    }
    return 0;
}

static int load_init(const char *fn, int drop) {
    FILE *f = fopen(fn, "r"); int n2; double sv, mdv;
    if (!f || fscanf(f, "%d %lf %lf", &n2, &sv, &mdv) != 3 || (n2 != N && n2 != N + 1)) { fprintf(stderr, "archivo inicial inválido\n"); return 0; }
    for (int i = 0; i < n2; i++) if (fscanf(f, "%lf %lf %lf", &X[i], &Y[i], &T[i]) != 3) return 0;
    fclose(f);
    if (n2 == N + 1) { X[drop] = X[N]; Y[drop] = Y[N]; T[drop] = T[N]; }
    S = sv; update_trig();
    return 1;
}

int main(int argc, char **argv) {
    if (argc == 2 && strcmp(argv[1], "gradtest") == 0) {
        N = 12; S = 5.0; rng_state = 12345;
        for (int i = 0; i < N; i++) { X[i] = 1 + 3 * urand(); Y[i] = 0.2 + 2.5 * urand(); T[i] = urand() * 3; }
        update_trig();
        static double g1[3 * MAXN], g2[3 * MAXN];
        gradient(g1); gradient_fd(g2);
        double md = 0, mx = 0;
        for (int k = 0; k < 3 * N; k++) { if (fabs(g1[k] - g2[k]) > md) md = fabs(g1[k] - g2[k]); if (fabs(g2[k]) > mx) mx = fabs(g2[k]); }
        printf("dif max analitico vs numerico = %.3e (escala %.3e), E=%.4f\n", md, mx, total_energy(NULL));
        return 0;
    }
    if (argc < 5) { fprintf(stderr, "uso: %s n segundos semilla s_inicial [inicial]\n", argv[0]); return 1; }
    N = atoi(argv[1]);
    budget = atof(argv[2]);
    seed = strtoull(argv[3], NULL, 10);
    double s0 = atof(argv[4]);
    const char *initfile = argc > 5 ? argv[5] : NULL;
    int fixed_mode = initfile && strcmp(initfile, "fixed") == 0;
    if (fixed_mode) { initfile = NULL; if (argc > 6) fixed_init = argv[6]; }
    rng_state = 0x9E3779B97F4A7C15ULL ^ (seed * 0xBF58476D1CE4E5B9ULL); for (int k = 0; k < 10; k++) urand();
    snprintf(out, sizeof out, "best_n%d_seed%llu.txt", N, seed);
    t0 = now();
    int restarts = 0, MBH_MAX = 60 + 4 * N;
    static double cx[MAXN], cy[MAXN], ct[MAXN];

    if (fixed_mode) {
        int tries = 0;
        while (now() - t0 < budget) { tries++; if (fixed_search(s0, 400 + 20 * N)) break; }
        printf("FIN n=%d seed=%llu mejor s = %.9f (intentos %d)\n", N, seed, best, tries);
        return 0;
    }
    while (now() - t0 < budget) {
        double cur;
        if (initfile && restarts == 0) {
            FILE *f = fopen(initfile, "r"); int n2; double sv, mdv;
            if (!f || fscanf(f, "%d %lf %lf", &n2, &sv, &mdv) != 3 || (n2 != N && n2 != N + 1)) { fprintf(stderr, "archivo inicial inválido\n"); return 1; }
            for (int i = 0; i < n2; i++) if (fscanf(f, "%lf %lf %lf", &X[i], &Y[i], &T[i]) != 3) return 1;
            fclose(f);
            if (n2 == N + 1) {
                int k = (int)(seed % (unsigned long long)n2);
                X[k] = X[N]; Y[k] = Y[N]; T[k] = T[N];
                printf("quitando el cuadrado %d de la configuración de %d\n", k, n2);
            }
            update_trig();
            cur = compress(sv, 1);
        } else {
            random_init(s0);
            cur = compress(s0, 1);
        }
        restarts++;
        if (cur < 0) continue;
        for (int i = 0; i < N; i++) { cx[i] = gx[i]; cy[i] = gy[i]; ct[i] = gt[i]; }
        // basin hopping
        int nofail = 0;
        while (nofail < MBH_MAX && now() - t0 < budget) {
            for (int i = 0; i < N; i++) { X[i] = cx[i]; Y[i] = cy[i]; T[i] = ct[i]; }
            S = cur; perturb(cur);
            double r = compress(cur * (1 - 1e-4), 0);
            if (r > 0 && r < cur - 1e-9) {
                cur = r; nofail = 0;
                for (int i = 0; i < N; i++) { cx[i] = gx[i]; cy[i] = gy[i]; ct[i] = gt[i]; }
            } else nofail++;
        }
    }
    printf("FIN n=%d seed=%llu mejor s = %.9f (reinicios %d)\n", N, seed, best, restarts);
    return 0;
}
