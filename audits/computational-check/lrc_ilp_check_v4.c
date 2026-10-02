/*
 * lrc_ilp_check_v4.c -- hunt mode: the probe-lemma filtered census.
 *
 * v4 = v3 + a "hunt" mode for pushing the exhaustive census to larger V.
 * solve(), the enumeration order, the record format, the tight/rung2/
 * interrung emission and the sample/spectrum machinery are UNTOUCHED
 * from v3; hunt mode adds:
 *
 *   hunt n V prefix S I tnum tden [nofilter]
 *
 *   Probe-lemma filter.  For any fixed rational probe t0 = a/p,
 *     gap(v) = sup_t min_j ||v_j t||  >=  min_j ||v_j (a/p)||.
 *   If min_j ||v_j (a/p)|| > theta (theta = tnum/tden), the vector's gap
 *   is CERTIFIED > theta and the full solve is skipped.  The filter is
 *   exact integer arithmetic; it can only skip vectors with gap > theta,
 *   so a census of all values <= theta is complete.  Sampled vectors
 *   (every 99733rd, as in v3) are always fully solved, preserving the
 *   S2 sample validation.
 *
 *   Extra outputs:
 *     <prefix>_lowgap.txt  every fully solved vector with gap <= theta
 *                          (v=(...) t=a/b gap=n/d): the exact low-gap
 *                          census at this radius;
 *     <prefix>_rungk.txt   k = 3..floor(V/n): vectors with gap exactly
 *                          k/(kn+1) (the rung-k value), same format as
 *                          the rung2 file;
 *     <prefix>_spectrum.txt  exact reduced histogram of ALL fully solved
 *                          gaps, plus a final line "filtered <count>"
 *                          (vectors certified gap > theta by a probe).
 *     <prefix>_summary.txt  total / filtered / fullsolved counters.
 *
 * With theta = 2/15 at n = 8 the filter keeps exact: tight (1/9), the
 * inter-rung window (1/9, 2/17), rung-2 (2/17), and every rung value
 * k/(8k+1) for k = 1..6 (all < 2/15), at a fraction of the full cost.
 *
 * Build: gcc -O2 -o lrc_ilp_check_v4 lrc_ilp_check_v4.c
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

/* ------------------------- shared state ------------------------- */
static FILE  *g_dump, *g_tight, *g_rung2, *g_interrung;
static int    g_n;
static long long g_total, g_feas, g_ntight, g_vfail, g_nrung2, g_ninter;
static long long g_r2_num, g_r2_den;      /* rung-2 target 2/(2n+1)  */

#define DENMAX 512                  /* max gap denominator: 2V, V <= 255 */
static long long g_tally[DENMAX + 1][DENMAX / 2 + 1];
static int    g_tally_on = 0;       /* spectrum histogram on/off        */
static long long g_sample_every = 0; /* 0 = no sample records           */
static FILE  *g_sample = NULL;

static int gcd_int(int a, int b)
{
    while (b) { int t = a % b; a = b; b = t; }
    return a;
}

/* tally the gap value num/den (unreduced) into the reduced histogram */
static void tally_gap(int num, int den)
{
    if (den < 1 || den > DENMAX) return;
    int g = gcd_int(num, den);
    if (g > 0) { num /= g; den /= g; }
    if (num > den / 2) return;       /* cannot happen: m <= b/2 in solve */
    g_tally[den][num]++;
}

/* C(m, r), exact in long long for m <= 255, r <= 8 */
static long long comb(int m, int r)
{
    if (r < 0 || m < r) return 0;
    long long acc = 1;
    for (int i = 0; i < r; i++)
        acc = acc * (m - i) / (i + 1);   /* exact at every step */
    return acc;
}

/* the r-th (0-based) vector in the lexicographic order of recurse() */
static void unrank_comb(long long r, int n, int V, int *v)
{
    int pos = 0, x = 1;
    while (pos < n) {
        if (x > V) { fprintf(stderr, "unrank error\n"); exit(1); }
        long long cnt = comb(V - x, n - 1 - pos);
        if (r < cnt) { v[pos++] = x; x++; }
        else { r -= cnt; x++; }
    }
}

/* advance v to the next vector of that same order; 0 if v was the last */
static int next_comb(int n, int V, int *v)
{
    int p = n - 1;
    while (p >= 0 && v[p] == V - (n - 1 - p)) p--;
    if (p < 0) return 0;
    v[p]++;
    for (int q = p + 1; q < n; q++) v[q] = v[q - 1] + 1;
    return 1;
}

/* ------------------------- hunt mode (v4) ------------------------- */

static long long g_filtered, g_fullsolved;
static FILE *g_lowgap = NULL;
static FILE *g_rungk[9];            /* k = 3..8: rung-k zoo files   */
static long long g_nrungk[9];
static int g_tnum, g_tden;          /* theta = tnum/tden            */
static int g_filter_on = 1;
static int g_K;                     /* emit rung k = 3..K           */

/*
 * The probe-lemma filter, v2: an early-exit scan over the solver's own
 * candidate times (peaks and crossings, the exact candidate set of
 * solve()).  If ANY candidate t = a/b gives min_j ||v_j t|| > theta,
 * then gap(v) >= that value > theta and the vector is certified out of
 * every hunted band (all values <= theta).  Typical vectors hit a
 * passing candidate within a handful of evaluations, so the scan costs
 * a few percent of a full solve; vectors with gap <= theta scan the
 * whole list (bounded by one full solve) and fall through to the exact
 * emission -- which is exactly the hunted population.
 */
static int filter_certifies(int n, const int *v)
{
    /* peaks: t = a/(2 v_i) */
    for (int i = 0; i < n; i++) {
        int b = 2 * v[i];
        for (int a = 1; a < b; a += 2) {
            int m = INT_MAX;
            for (int j = 0; j < n; j++) {
                int r = (int)((long long)a * v[j] % b);
                int d = r < b - r ? r : b - r;
                if (d < m) m = d;
            }
            if ((long long)m * g_tden > (long long)g_tnum * b)
                return 1;
        }
    }
    /* crossings: t = a/(v_i+v_j), a/(v_j-v_i) */
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            int b = v[i] + v[j];
            for (int a = 0; a <= b; a++) {
                int m = INT_MAX;
                for (int q = 0; q < n; q++) {
                    int r = (int)((long long)a * v[q] % b);
                    int d = r < b - r ? r : b - r;
                    if (d < m) m = d;
                }
                if ((long long)m * g_tden > (long long)g_tnum * b)
                    return 1;
            }
            b = v[j] - v[i];
            for (int a = 0; a <= b; a++) {
                int m = INT_MAX;
                for (int q = 0; q < n; q++) {
                    int r = (int)((long long)a * v[q] % b);
                    int d = r < b - r ? r : b - r;
                    if (d < m) m = d;
                }
                if ((long long)m * g_tden > (long long)g_tnum * b)
                    return 1;
            }
        }
    return 0;
}

/* full emission for one vector (v3 emit() + v4 extras); n = g_n. */
static void emit_full(const int *v)
{
    int k[NVMAX], ta, tb, num, den;
    g_fullsolved++;

    solve(g_n, v, k, &ta, &tb, &num, &den);
    int feas = ((long long)num * (g_n + 1) >= (long long)den);
    int is_r2  = ((long long)num * (g_r2_den) == 2 * (long long)den);
    int is_ir  = ((long long)num * (g_n + 1) > (long long)den) &&
                 ((long long)num * (g_r2_den) < 2 * (long long)den);
    long long gap_le_theta =
        (long long)num * g_tden <= (long long)g_tnum * den;

    if (g_tally_on) tally_gap(num, den);
    if (g_sample && g_sample_every && g_total % g_sample_every == 0) {
        fprintf(g_sample, "v=(");
        for (int i = 0; i < g_n; i++)
            fprintf(g_sample, "%d%s", v[i], i + 1 < g_n ? "," : "");
        fprintf(g_sample, ") k=(");
        for (int i = 0; i < g_n; i++)
            fprintf(g_sample, "%d%s", k[i], i + 1 < g_n ? "," : "");
        fprintf(g_sample, ") t=%d/%d gap=%d/%d feasible=%d\n",
                ta, tb, num, den, feas);
    }
    if (gap_le_theta && g_lowgap) {
        fprintf(g_lowgap, "v=(");
        for (int i = 0; i < g_n; i++)
            fprintf(g_lowgap, "%d%s", v[i], i + 1 < g_n ? "," : "");
        fprintf(g_lowgap, ") t=%d/%d gap=%d/%d\n", ta, tb, num, den);
    }
    for (int kk = 3; kk <= g_K; kk++) {
        long long kn = (long long)kk * g_n + 1;
        if ((long long)num * kn == (long long)kk * den && g_rungk[kk]) {
            g_nrungk[kk]++;
            fprintf(g_rungk[kk], "v=(");
            for (int i = 0; i < g_n; i++)
                fprintf(g_rungk[kk], "%d%s", v[i], i + 1 < g_n ? "," : "");
            fprintf(g_rungk[kk], ") t=%d/%d gap=%d/%d bind=(", ta, tb, num, den);
            int m = INT_MAX;
            for (int i = 0; i < g_n; i++) {
                int rr = (int)((long long)ta * v[i] % tb);
                int dd = rr < tb - rr ? rr : tb - rr;
                if (dd < m) m = dd;
            }
            int first = 1;
            for (int i = 0; i < g_n; i++) {
                int rr = (int)((long long)ta * v[i] % tb);
                int dd = rr < tb - rr ? rr : tb - rr;
                if (dd == m) {
                    fprintf(g_rungk[kk], "%s%d", first ? "" : ",", v[i]);
                    first = 0;
                }
            }
            fprintf(g_rungk[kk], ")\n");
        }
    }

    g_total++;
    if (is_r2) {
        g_nrung2++;
        fprintf(g_rung2, "v=(");
        for (int i = 0; i < g_n; i++)
            fprintf(g_rung2, "%d%s", v[i], i+1<g_n?",":"");
        fprintf(g_rung2, ") t=%d/%d gap=%d/%d bind=(", ta, tb, num, den);
        int m = INT_MAX;
        for (int i = 0; i < g_n; i++) {
            int rr = (int)((long long)ta * v[i] % tb);
            int dd = rr < tb - rr ? rr : tb - rr;
            if (dd < m) m = dd;
        }
        int first = 1;
        for (int i = 0; i < g_n; i++) {
            int rr = (int)((long long)ta * v[i] % tb);
            int dd = rr < tb - rr ? rr : tb - rr;
            if (dd == m) {
                fprintf(g_rung2, "%s%d", first ? "" : ",", v[i]);
                first = 0;
            }
        }
        fprintf(g_rung2, ")\n");
    }
    if (is_ir) {
        g_ninter++;
        fprintf(g_interrung, "v=(");
        for (int i = 0; i < g_n; i++)
            fprintf(g_interrung, "%d%s", v[i], i+1<g_n?",":"");
        fprintf(g_interrung, ") t=%d/%d gap=%d/%d\n", ta, tb, num, den);
    }
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
            for (int i = 0; i < g_n; i++)
                fprintf(g_tight, "%d%s", v[i], i+1<g_n?",":"");
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

static void run_hunt_mode(int n, int V, const char *prefix,
                          int shards, int idx, int dump_on,
                          int tnum, int tden, int filter_on)
{
    char path[1024];
    int v[NVMAX];
    long long total, base, rem, lo, hi, cnt, done;

    if (n < 1 || n > NVMAX) { fprintf(stderr, "bad n\n"); exit(1); }
    if (V < n || V > 255)   { fprintf(stderr, "bad V\n"); exit(1); }
    if (shards < 1 || idx < 0 || idx >= shards) {
        fprintf(stderr, "bad shard index\n"); exit(1);
    }

    g_n = n; g_total = g_feas = g_ntight = g_vfail = 0;
    g_nrung2 = g_ninter = 0;
    g_r2_num = 2; g_r2_den = 2 * n + 1;
    g_tally_on = 1; g_sample_every = 99733;
    g_filtered = g_fullsolved = 0;
    g_tnum = tnum; g_tden = tden; g_filter_on = filter_on;
    g_K = V / n;
    for (int kk = 3; kk <= 8; kk++) { g_rungk[kk] = NULL; g_nrungk[kk] = 0; }

    total = comb(V, n);
    base = total / shards;
    rem = total % shards;
    lo = (long long)idx * base + (idx < rem ? idx : rem);
    hi = (long long)(idx + 1) * base + (idx + 1 < rem ? idx + 1 : rem);
    cnt = hi - lo;

    snprintf(path, sizeof path, "%s_dump.bin", prefix);
    g_dump = dump_on ? fopen(path, "wb") : NULL;
    snprintf(path, sizeof path, "%s_tight.txt", prefix);
    g_tight = fopen(path, "w");
    snprintf(path, sizeof path, "%s_rung2.txt", prefix);
    g_rung2 = fopen(path, "w");
    snprintf(path, sizeof path, "%s_interrung.txt", prefix);
    g_interrung = fopen(path, "w");
    snprintf(path, sizeof path, "%s_sample.txt", prefix);
    g_sample = fopen(path, "w");
    snprintf(path, sizeof path, "%s_lowgap.txt", prefix);
    g_lowgap = fopen(path, "w");
    for (int kk = 3; kk <= g_K && kk <= 8; kk++) {
        snprintf(path, sizeof path, "%s_rung%d.txt", prefix, kk);
        g_rungk[kk] = fopen(path, "w");
    }
    if ((dump_on && !g_dump) || !g_tight || !g_rung2 || !g_interrung
        || !g_sample || !g_lowgap) {
        perror("fopen"); exit(1);
    }

    clock_t t0 = clock();
    unrank_comb(lo, n, V, v);
    for (done = 0; done < cnt; done++) {
        if ((done & 0x1FFFFF) == 0)
            fprintf(stderr, "shard %d/%d: %lld / %lld\n",
                    idx, shards, done, cnt);
        if (g_filter_on
            && !(g_sample_every && g_total % g_sample_every == 0)
            && filter_certifies(n, v)) {
            g_filtered++;
            g_total++;
        } else {
            emit_full(v);
        }
        if (done + 1 < cnt && !next_comb(n, V, v)) {
            fprintf(stderr, "shard enumeration ended early\n"); exit(1);
        }
    }
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;

    if (g_dump) fclose(g_dump);
    fclose(g_tight); fclose(g_rung2); fclose(g_interrung);
    fclose(g_sample); fclose(g_lowgap);
    for (int kk = 3; kk <= 8; kk++)
        if (g_rungk[kk]) { fclose(g_rungk[kk]); }

    snprintf(path, sizeof path, "%s_spectrum.txt", prefix);
    FILE *fsp = fopen(path, "w");
    if (!fsp) { perror("fopen"); exit(1); }
    for (int den = 1; den <= DENMAX; den++)
        for (int num = 0; num <= den / 2; num++)
            if (g_tally[den][num])
                fprintf(fsp, "%d %d %lld\n", num, den, g_tally[den][num]);
    fprintf(fsp, "filtered %lld\n", g_filtered);
    fclose(fsp);

    snprintf(path, sizeof path, "%s_summary.txt", prefix);
    FILE *fs = fopen(path, "w");
    fprintf(fs, "n=%d V=%d shards=%d idx=%d lo=%lld hi=%lld count=%lld "
                "filtered=%lld fullsolved=%lld feasible=%lld "
                "infeasible=%lld tight=%lld rung2=%lld interrung=%lld ",
            n, V, shards, idx, lo, hi, cnt, g_filtered, g_fullsolved,
            g_feas, g_fullsolved - g_feas, g_ntight, g_nrung2, g_ninter);
    for (int kk = 3; kk <= g_K && kk <= 8; kk++)
        fprintf(fs, "rung%d=%lld ", kk, g_nrungk[kk]);
    fprintf(fs, "verify_failures=%lld theta=%d/%d seconds=%.2f\n",
            g_vfail, g_tnum, g_tden, secs);
    fclose(fs);
    printf("n=%d V=%d shards=%d idx=%d count=%lld filtered=%lld "
           "fullsolved=%lld tight=%lld rung2=%lld interrung=%lld ",
           n, V, shards, idx, cnt, g_filtered, g_fullsolved,
           g_ntight, g_nrung2, g_ninter);
    for (int kk = 3; kk <= g_K && kk <= 8; kk++)
        printf("rung%d=%lld ", kk, g_nrungk[kk]);
    printf("seconds=%.2f\n", secs);
}

int main(int argc, char **argv)
{
    if (argc >= 2 && !strcmp(argv[1], "run") && argc == 5) {
        /* v3 run mode unchanged (not used by the hunt) */
        fprintf(stderr, "run mode: use lrc_ilp_check_v3\n");
        return 1;
    } else if (argc >= 2 && !strcmp(argv[1], "hunt")) {
        /* hunt n V prefix S I tnum tden [nodump] [nofilter] */
        int dump_on = 1, filter_on = 1;
        if (argc < 9) {
            fprintf(stderr,
                    "usage: %s hunt <n> <V> <prefix> <shards> <idx> "
                    "<tnum> <tden> [nodump] [nofilter]\n", argv[0]);
            return 1;
        }
        for (int i = 9; i < argc; i++) {
            if (!strcmp(argv[i], "nodump")) dump_on = 0;
            if (!strcmp(argv[i], "nofilter")) filter_on = 0;
        }
        run_hunt_mode(atoi(argv[2]), atoi(argv[3]), argv[4],
                      atoi(argv[5]), atoi(argv[6]), dump_on,
                      atoi(argv[7]), atoi(argv[8]), filter_on);
        return 0;
    } else if (argc >= 2 && !strcmp(argv[1], "shard") &&
               (argc == 7 || (argc == 8 && !strcmp(argv[7], "nodump")))) {
        fprintf(stderr, "shard mode: use lrc_ilp_check_v3\n");
        return 1;
    } else if (argc >= 2 && !strcmp(argv[1], "brute") && argc == 3) {
        fprintf(stderr, "brute mode: use lrc_ilp_check_v3\n");
        return 1;
    } else {
        fprintf(stderr,
                "usage: %s hunt <n> <V> <prefix> <shards> <idx> "
                "<tnum> <tden> [nodump] [nofilter]\n", argv[0]);
        return 1;
    }
}
