/*
 * lrc_ilp_check.c
 *
 * Exhaustive integer-feasibility check for the Lonely Runner polytope:
 *
 *   v_1 < ... < v_n distinct positive integers, delta = 1/(n+1).
 *   P(v) = { k in R^n : 0 <= k_i <= v_i - 1,
 *            (n+1)(v_j k_i - v_i k_j) <= n v_i - v_j   for all i != j }.
 *
 *   Conjecture (equivalent form): P(v) cap Z^n != empty for every v.
 *
 * Method (exact, rational arithmetic simulated with integers):
 *   gap(v) = max_{t in [0,1]} min_i || t v_i ||_T .
 *   The function g(t) = min_i ||t v_i|| is piecewise linear and continuous;
 *   its maximum over [0,1] is attained at a breakpoint of some ||t v_i||
 *   (t = odd/(2 v_i), the "peaks") or at a crossing ||t v_i|| = ||t v_j||
 *   (t(v_i+v_j) in Z or t(v_i-v_j) in Z), or at t in {0,1}.  Valleys
 *   (t = m/v_i) give g = 0 and can never host a positive maximum.
 *
 *   gap(v) >= 1/(n+1)  <=>  P(v) contains an integer point.  A witness is
 *   k_i = floor(t* v_i) where t* is an argmax; it is verified against all
 *   polytope constraints (integer arithmetic) for every single vector.
 *
 * Modes:
 *   run   n V prefix   : enumerate all v with v_n <= V; write
 *                        <prefix>_dump.bin (32-byte records),
 *                        <prefix>_summary.txt; tight list to stdout path
 *                        <prefix>_tight.txt (gap == 1/(n+1) exactly).
 *   brute n            : read vectors "v1 ... vn" per line from stdin and
 *                        enumerate the ENTIRE k-box [0,v_i-1]^n, checking
 *                        every integer point against all constraints
 *                        (literal ILP enumeration).  Prints count of
 *                        feasible integer points and the first one found.
 *
 * Build: gcc -O2 -o lrc_ilp_check lrc_ilp_check.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <limits.h>

#define NVMAX   8
#define CANDMAX 8192

#pragma pack(push, 1)
typedef struct {
    uint8_t  v[NVMAX];   /* speeds, zero-padded                       */
    uint8_t  k[NVMAX];   /* witness integer point, zero-padded        */
    uint16_t ta, tb;     /* argmax t = ta/tb                          */
    uint16_t num, den;   /* gap = num/den (unreduced)                 */
    uint8_t  feasible;   /* 1 iff gap*(n+1) >= 1                      */
    uint8_t  pad[7];
} Rec;
#pragma pack(pop)

/* Verify an integer point k against ALL polytope constraints for speeds v. */
static int verify_point(int n, const int *v, const int *k)
{
    for (int i = 0; i < n; i++)
        if (k[i] < 0 || k[i] > v[i] - 1)
            return 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            long long lhs = (long long)(n + 1) *
                ((long long)v[j] * k[i] - (long long)v[i] * k[j]);
            long long rhs = (long long)n * v[i] - v[j];
            if (lhs > rhs)
                return 0;
        }
    }
    return 1;
}

/* Exact gap computation; outputs witness k, argmax t = ta/tb, gap = num/den. */
static void solve(int n, const int *v, int *k,
                  int *ta, int *tb, int *num, int *den)
{
    static int ca[CANDMAX], cb[CANDMAX];
    int nc = 0;

    /* peaks of ||t v_i||: t = odd/(2 v_i) */
    for (int i = 0; i < n; i++)
        for (int a = 1; a < 2 * v[i]; a += 2) {
            if (nc >= CANDMAX) { fprintf(stderr, "cand overflow\n"); exit(1); }
            ca[nc] = a; cb[nc] = 2 * v[i]; nc++;
        }
    /* crossings ||t v_i|| = ||t v_j||: t (v_i+v_j) in Z or t (v_j-v_i) in Z */
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            int s = v[i] + v[j];
            for (int a = 0; a <= s; a++) {
                if (nc >= CANDMAX) { fprintf(stderr, "cand overflow\n"); exit(1); }
                ca[nc] = a; cb[nc] = s; nc++;
            }
            int d = v[j] - v[i];
            for (int a = 0; a <= d; a++) {
                if (nc >= CANDMAX) { fprintf(stderr, "cand overflow\n"); exit(1); }
                ca[nc] = a; cb[nc] = d; nc++;
            }
        }

    int bn = 0, bd = 1, ba = 0, bb = 1;   /* best value 0 at t = 0 */
    for (int c = 0; c < nc; c++) {
        int a = ca[c], b = cb[c];
        int m = INT_MAX;
        for (int i = 0; i < n; i++) {
            int r = (int)((long long)a * v[i] % b);
            int d = r < b - r ? r : b - r;
            if (d < m) m = d;
        }
        /* value m/b vs best bn/bd */
        if ((long long)m * bd > (long long)bn * b) {
            bn = m; bd = b; ba = a; bb = b;
        }
    }
    for (int i = 0; i < n; i++)
        k[i] = (int)((long long)ba * v[i] / bb);
    *ta = ba; *tb = bb; *num = bn; *den = bd;
}

/* ------------------------- run mode ------------------------- */
static FILE  *g_dump, *g_tight;
static int    g_n;
static long long g_total, g_feas, g_ntight, g_vfail;

static void emit(const int *v)
{
    int k[NVMAX], ta, tb, num, den;
    Rec r;
    memset(&r, 0, sizeof r);

    solve(g_n, v, k, &ta, &tb, &num, &den);
    int feas = ((long long)num * (g_n + 1) >= (long long)den);

    for (int i = 0; i < g_n; i++) { r.v[i] = (uint8_t)v[i]; r.k[i] = (uint8_t)k[i]; }
    r.ta = (uint16_t)ta; r.tb = (uint16_t)tb;
    r.num = (uint16_t)num; r.den = (uint16_t)den;
    r.feasible = (uint8_t)feas;
    fwrite(&r, sizeof r, 1, g_dump);

    g_total++;
    if (feas) {
        g_feas++;
        if (!verify_point(g_n, v, k)) {
            g_vfail++;
            fprintf(stderr, "VERIFY FAILURE v=(");
            for (int i = 0; i < g_n; i++) fprintf(stderr, "%d%s", v[i], i+1<g_n?",":"");
            fprintf(stderr, ") k=(");
            for (int i = 0; i < g_n; i++) fprintf(stderr, "%d%s", k[i], i+1<g_n?",":"");
            fprintf(stderr, ")\n");
        }
        if ((long long)num * (g_n + 1) == (long long)den) {   /* tight case */
            g_ntight++;
            fprintf(g_tight, "v=(");
            for (int i = 0; i < g_n; i++) fprintf(g_tight, "%d%s", v[i], i+1<g_n?",":"");
            fprintf(g_tight, ") k=(");
            for (int i = 0; i < g_n; i++) fprintf(g_tight, "%d%s", k[i], i+1<g_n?",":"");
            fprintf(g_tight, ") t=%d/%d gap=%d/%d\n", ta, tb, num, den);
        }
    } else {
        fprintf(stderr, "INFEASIBLE v=(");
        for (int i = 0; i < g_n; i++) fprintf(stderr, "%d%s", v[i], i+1<g_n?",":"");
        fprintf(stderr, ") gap=%d/%d\n", num, den);
    }
}

static void recurse(int pos, int start, int V, int *v)
{
    if (pos == g_n) { emit(v); return; }
    for (int x = start; x <= V - (g_n - 1 - pos); x++) {
        v[pos] = x;
        recurse(pos + 1, x + 1, V, v);
    }
}

static void run_mode(int n, int V, const char *prefix)
{
    char path[1024];
    int v[NVMAX];
    g_n = n; g_total = g_feas = g_ntight = g_vfail = 0;

    snprintf(path, sizeof path, "%s_dump.bin", prefix);
    g_dump = fopen(path, "wb");
    snprintf(path, sizeof path, "%s_tight.txt", prefix);
    g_tight = fopen(path, "w");
    if (!g_dump || !g_tight) { perror("fopen"); exit(1); }

    clock_t t0 = clock();
    recurse(0, 1, V, v);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;

    fclose(g_dump); fclose(g_tight);

    snprintf(path, sizeof path, "%s_summary.txt", prefix);
    FILE *fs = fopen(path, "w");
    fprintf(fs, "n=%d V=%d total=%lld feasible=%lld infeasible=%lld "
                "tight=%lld verify_failures=%lld seconds=%.2f\n",
            n, V, g_total, g_feas, g_total - g_feas, g_ntight, g_vfail, secs);
    fclose(fs);
    printf("n=%d V=%d total=%lld feasible=%lld infeasible=%lld tight=%lld "
           "verify_failures=%lld seconds=%.2f\n",
           n, V, g_total, g_feas, g_total - g_feas, g_ntight, g_vfail, secs);
}

/* ------------------------- brute mode ------------------------- */
static void brute_mode(int n)
{
    int v[NVMAX], k[NVMAX], first[NVMAX];
    long long nvec = 0, nsol_total = 0, nvec_with_sol = 0;
    while (1) {
        int ok = 1;
        for (int i = 0; i < n; i++)
            if (scanf("%d", &v[i]) != 1) { ok = 0; break; }
        if (!ok) break;
        nvec++;
        long long nsol = 0; int found = 0;
        for (int i = 0; i < n; i++) k[i] = 0;
        for (;;) {
            if (verify_point(n, v, k)) {
                if (!found) { memcpy(first, k, sizeof(int) * n); found = 1; }
                nsol++;
            }
            int i = 0;
            while (i < n && k[i] == v[i] - 1) k[i++] = 0;
            if (i == n) break;
            k[i]++;
        }
        nsol_total += nsol;
        if (found) nvec_with_sol++;
        printf("v=(");
        for (int i = 0; i < n; i++) printf("%d%s", v[i], i+1<n?",":"");
        printf(") feasible_integer_points=%lld first_k=(", nsol);
        if (found)
            for (int i = 0; i < n; i++) printf("%d%s", first[i], i+1<n?",":"");
        else printf("NONE");
        printf(")\n");
    }
    fprintf(stderr, "brute: vectors=%lld with_solution=%lld total_solutions=%lld\n",
            nvec, nvec_with_sol, nsol_total);
}

int main(int argc, char **argv)
{
    if (argc >= 2 && !strcmp(argv[1], "run") && argc == 5) {
        run_mode(atoi(argv[2]), atoi(argv[3]), argv[4]);
    } else if (argc >= 2 && !strcmp(argv[1], "brute") && argc == 3) {
        brute_mode(atoi(argv[2]));
    } else {
        fprintf(stderr,
                "usage: %s run <n> <V> <prefix>\n"
                "       %s brute <n>   (vectors on stdin)\n", argv[0], argv[0]);
        return 1;
    }
    return 0;
}
