/*
 * lrc_xs_batch.c -- batch X(S) brute force for the filler-characterization
 * verification (lrc_xs_theory.py).
 *
 * usage: lrc_xs_batch <gnum> <gden> <xmax>
 * stdin: one core per line:  "k v1 v2 ... vk"   (k = n-1 speeds, sorted)
 * stdout per core:
 *   "v1,v2,...,vk | coregap=num/den | X=x1,x2,..."
 * where X = { x in [1,xmax], x not in core : gap(core + x) = gnum/gden }.
 *
 * The gap computation is solve() from lrc_ilp_check_v3.c, VERBATIM
 * (same candidate set, same integer arithmetic, same reduction-free
 * comparison), so the brute-force side is the sweep solver's own
 * arithmetic; the theory side (lrc_xs_theory.py) is independent Python
 * Fraction code.
 *
 * Build: gcc -O2 -o lrc_xs_batch lrc_xs_batch.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>

#define NVMAX   16
#define CANDMAX 8192

/* ---------------- solve(), copied verbatim from v3 ---------------- */
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
/* ------------------------------------------------------------------- */

int main(int argc, char **argv)
{
    if (argc != 4) {
        fprintf(stderr, "usage: %s <gnum> <gden> <xmax>\n", argv[0]);
        return 1;
    }
    long long gnum = atoll(argv[1]), gden = atoll(argv[2]);
    int xmax = atoi(argv[3]);

    char line[4096];
    long long ncores = 0;
    while (fgets(line, sizeof line, stdin)) {
        int v[NVMAX], w[NVMAX + 1], k[NVMAX + 1];
        int kargs = 0, off = 0, read;
        int nc;
        if (sscanf(line + off, "%d%n", &nc, &read) != 1) continue;
        off += read;
        for (int i = 0; i < nc; i++) {
            if (sscanf(line + off, "%d%n", &v[i], &read) != 1) {
                fprintf(stderr, "bad line: %s", line); return 1;
            }
            off += read;
        }
        int ta, tb, num, den;
        solve(nc, v, k, &ta, &tb, &num, &den);
        long long cnum = num, cden = den;

        /* print core and its gap */
        for (int i = 0; i < nc; i++)
            printf("%s%d", i ? "," : "", v[i]);
        printf(" | %lld/%lld | ", cnum, cden);

        int first = 1;
        for (int x = 1; x <= xmax; x++) {
            int incore = 0;
            for (int i = 0; i < nc; i++) if (v[i] == x) { incore = 1; break; }
            if (incore) continue;
            for (int i = 0; i < nc; i++) w[i] = v[i];
            w[nc] = x;
            /* keep sorted order not required by solve(); keep as-is */
            solve(nc + 1, w, k, &ta, &tb, &num, &den);
            if ((long long)num * gden == gnum * (long long)den) {
                printf("%s%d", first ? "" : ",", x);
                first = 0;
            }
        }
        printf("\n");
        ncores++;
    }
    fprintf(stderr, "cores=%lld\n", ncores);
    return 0;
}
