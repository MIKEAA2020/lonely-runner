/*
 * pool_sat.c -- saturation checker for the frontier-cell pool law.  v2
 *
 * QUESTION (reviewer, 2026-10-06): is the pool identity
 *     n_distinct_families = (D-1)(D-2)...(D-(m-1)),  D = (p-1)/2
 * (measured EXACT at T=8 (1K,5U) p=17: 840 = 7*6*5*4) a theorem or a
 * p=17 coincidence?
 *
 * For each ORDERED DISTINCT difference tuple (d2..dm) with values in
 * [2,D], pairwise distinct, decide ADMISSIBILITY with the committed
 * census_m5 convention (non-bystander):
 *   admissible  <=>  exists ordered size tuple (s1..sm) in {on,off}^m
 *   with sum >= p, and starts a2..a_{m-1} (a1=0, d1=1) with arcs
 *   1..m-1 placed, remainder R = Z_p \ union with 1 <= |R| <= s_m,
 *   and arc m (difference d_m, size s_m) containing R (arc-fit test).
 * R empty (arcs 1..m-1 already cover) is counted separately as
 * BYSTANDER-admissible (census_m5 excludes those).
 *
 * Pruning mirrors census_m5 exactly, rewritten popcount-lean:
 *   inter = |union_prev & arc| <= ov, and
 *   |union_prev| + s_i - inter >= p - s_{i+1..m}   (uncovered budget).
 *
 * Usage:  pool_sat T p m [cap]     (frontier cell (1K,(m-1)U), m = T-3)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define MAXM 8

typedef struct {
    int p, T, m, k, eps, on, off, D;
    uint64_t full;
} Ctx;

static Ctx C;
static uint64_t base_on[64], base_off[64];
static int invd[64];

static int d[MAXM + 1], s[MAXM + 1];
static uint64_t bm[MAXM + 1];
static uint64_t rotm[64];      /* rotations of arc (m-1)'s base mask */
static int ov;
static long long steps, step_cap;
static int found, bystander, capped;

static inline int popc(uint64_t x) { return __builtin_popcountll(x); }

static inline uint64_t rot64(uint64_t x, int a, int p, uint64_t full) {
    if (a == 0) return x & full;
    return (((x << a) | (x >> (p - a))) & full);
}

static uint64_t ap_base(int d, int sz, int p) {
    uint64_t m = 0;
    long long x = 0;
    for (int j = 0; j < sz; j++) { m |= 1ULL << x; x += d; if (x >= p) x -= p; }
    return m;
}

/* arc-fit test (mirrors fit_placements): can the AP with difference dm,
 * size sm contain all residues of R (|R| = nr in [1, sm])? */
static int fit_test(uint64_t R, int dm, int sm) {
    int v[64], nv = 0;
    int iv = invd[dm];
    uint64_t r = R;
    while (r) {
        int b = __builtin_ctzll(r);
        v[nv++] = (int)((long long)iv * b % C.p);
        r &= r - 1;
    }
    int arc_len;
    if (nv == 1) arc_len = 1;
    else {
        for (int i = 1; i < nv; i++) { int x = v[i], j = i - 1;
            while (j >= 0 && v[j] > x) { v[j + 1] = v[j]; j--; } v[j + 1] = x; }
        int gmax = 0;
        for (int i = 0; i < nv; i++) {
            int g = v[(i + 1) % nv] - v[i];
            if (g < 0) g += C.p;
            if (g > gmax) gmax = g;
        }
        arc_len = C.p - gmax + 1;
    }
    return arc_len <= sm;
}

/* levels 2..m-2 recursive; level m-1 flattened; level m = fit test */
static void rec(int i, uint64_t u, int pu) {
    if (found || capped) return;
    if (i == C.m - 1) {                 /* flattened deepest loop */
        int s4 = s[i], s5 = s[C.m];
        int dm = d[C.m];
        for (int a = 0; a < C.p && !found && !capped; a++) {
            if (step_cap > 0 && ++steps > step_cap) { capped = 1; return; }
            uint64_t mk = rotm[a];
            int inter = popc(u & mk);
            if (inter > ov) continue;
            int cov = pu + s4 - inter;            /* |u | mk| */
            if (cov < C.p - s5) continue;         /* uncovered > s5 */
            if (cov >= C.p) { bystander = 1; continue; }
            /* 1 <= nr <= s5 here; nr = p - cov */
            if (fit_test(C.full & ~(u | mk), dm, s5)) { found = 1; return; }
        }
        return;
    }
    uint64_t b = bm[i];
    int rem = 0;
    for (int j = i + 1; j <= C.m; j++) rem += s[j];
    for (int a = 0; a < C.p && !found && !capped; a++) {
        if (step_cap > 0 && ++steps > step_cap) { capped = 1; return; }
        uint64_t mk = rot64(b, a, C.p, C.full);
        int inter = popc(u & mk);
        if (inter > ov) continue;
        int cov = pu + s[i] - inter;
        if (C.p - cov > rem) continue;            /* uncovered > remaining */
        rec(i + 1, u | mk, cov);
    }
}

typedef struct { int s[MAXM + 1]; int sum; } SzT;
static SzT szts[1 << MAXM];
static int n_szt;

static void build_size_tuples(void) {
    n_szt = 0;
    for (int code = 0; code < (1 << C.m); code++) {
        SzT t; t.sum = 0;
        for (int i = 1; i <= C.m; i++) {
            t.s[i] = ((code >> (i - 1)) & 1) ? C.off : C.on;
            t.sum += t.s[i];
        }
        if (t.sum >= C.p) szts[n_szt++] = t;
    }
    for (int i = 1; i < n_szt; i++) { SzT x = szts[i]; int j = i - 1;
        while (j >= 0 && szts[j].sum < x.sum) { szts[j + 1] = szts[j]; j--; }
        szts[j + 1] = x; }
}

static long long n_tuples, n_adm, n_byst, n_unres;
static FILE *dumpf = NULL;
static long long skip_count = 0;   /* resume: tuples already done  */
static long long done_count = 0;
static long long visited = 0;      /* all tuples visited (pre-stride) */
static long long sample_stride = 0; /* 0 = run all; else run every stride-th */

static void do_tuple(void) {
    visited++;
    if (sample_stride > 0 && (visited % sample_stride) != 0) return;
    if (done_count < skip_count) { done_count++; return; }
    done_count++;
    found = 0; bystander = 0; capped = 0; steps = 0;
    for (int t = 0; t < n_szt && !found && !capped; t++) {
        for (int i = 1; i <= C.m; i++) s[i] = szts[t].s[i];
        ov = szts[t].sum - C.p;
        for (int i = 2; i <= C.m; i++)
            bm[i] = (s[i] == C.off) ? base_off[d[i]] : base_on[d[i]];
        for (int a = 0; a < C.p; a++) rotm[a] = rot64(bm[C.m - 1], a, C.p, C.full);
        uint64_t m1 = (1ULL << s[1]) - 1ULL;
        rec(2, m1 & C.full, s[1]);
    }
    n_tuples++;
    if (capped)         n_unres++;
    else if (found)     n_adm++;
    else if (bystander) n_byst++;
    if (dumpf) {
        fprintf(dumpf, "%d", d[2]);
        for (int i = 3; i <= C.m; i++) fprintf(dumpf, ",%d", d[i]);
        fprintf(dumpf, "\t%d\t%lld\n",
                capped ? -1 : (found ? 1 : (bystander ? 2 : 0)), steps);
        fflush(dumpf);
    }
}

static void gen(int pos, uint32_t used) {
    if (pos == C.m + 1) { do_tuple(); return; }
    for (int v = 2; v <= C.D; v++) {
        if (used & (1u << v)) continue;
        d[pos] = v;
        gen(pos + 1, used | (1u << v));
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s T p m [cap]\n", argv[0]); return 1; }
    C.T = atoi(argv[1]); C.p = atoi(argv[2]); C.m = atoi(argv[3]);
    long long cap = (argc > 4) ? atoll(argv[4]) : 0;
    int p = C.p, T = C.T;
    C.eps = p % T;
    C.k = (p - C.eps) / T;
    int d0 = (C.eps * C.eps - 1) / T;
    int M = C.k * C.eps + d0;
    if (2 * M >= p - 1) { C.on = 2 * C.k + 1; C.off = 2 * C.k + 2; }
    else                { C.on = 2 * C.k;     C.off = 2 * C.k + 1; }
    C.D = (p - 1) / 2;
    C.full = (1ULL << p) - 1ULL;
    step_cap = cap;
    if (argc > 6) sample_stride = atoll(argv[6]);
    if (C.m > MAXM || C.m < 3 || C.D > 31 || p > 63) { fprintf(stderr, "out of range\n"); return 1; }

    for (int dd = 1; dd <= C.D; dd++) {
        base_on[dd] = ap_base(dd, C.on, p);
        base_off[dd] = ap_base(dd, C.off, p);
        invd[dd] = -1;
        for (int x = 1; x < p; x++) if ((long long)dd * x % p == 1) { invd[dd] = x; break; }
    }
    build_size_tuples();
    if (argc > 5) {
        /* resume if the dump already exists; recount its decisions */
        FILE *probe = fopen(argv[5], "r");
        if (probe) {
            char buf[256];
            while (fgets(buf, sizeof buf, probe)) {
                if (buf[0] < '0' || buf[0] > '9') continue;
                skip_count++; n_tuples++;
                char *tab = strchr(buf, '\t');
                if (!tab) continue;
                int dec = atoi(tab + 1);
                if (dec == 1) n_adm++;
                else if (dec == 2) n_byst++;
                else if (dec == -1) n_unres++;
            }
            fclose(probe);
            dumpf = fopen(argv[5], "a");
        } else {
            dumpf = fopen(argv[5], "w");
        }
    }

    clock_t t0 = clock();
    gen(2, 0u);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    if (dumpf) { fclose(dumpf); dumpf = NULL; }
    if (skip_count)
        printf("(resumed after %lld tuples) ", skip_count);

    long long ff = 1;
    for (int j = 1; j <= C.m - 1; j++) ff *= (C.D - j);
    if (sample_stride > 0)
        printf("SAMPLE stride=%lld of %lld: ", sample_stride, ff);
    printf("{\"T\":%d,\"p\":%d,\"m\":%d,\"k\":%d,\"eps\":%d,\"on\":%d,\"off\":%d,"
           "\"D\":%d,\"n_size_tuples\":%d,\"step_cap\":%lld,"
           "\"n_tuples\":%lld,\"falling_factorial\":%lld,"
           "\"admissible\":%lld,\"bystander_only\":%lld,\"unresolved\":%lld,"
           "\"saturated\":%s,\"seconds\":%.2f}\n",
           T, p, C.m, C.k, C.eps, C.on, C.off, C.D, n_szt, cap,
           n_tuples, ff, n_adm, n_byst, n_unres,
           (n_unres == 0 && n_adm == ff) ? "true" :
           (n_unres == 0 && n_adm < ff ? "false" : "undetermined"),
           secs);
    fflush(stdout);
    return 0;
}
