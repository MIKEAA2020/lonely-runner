/* lrc_nsk_census.c -- the no-silent-kill (NSK) census.
 *
 * Background (audits: "rung family and the filler theorem.txt",
 * "tight-set decomposition via X(S).txt", "silent killers and the
 * induction step.txt").  General Filler Theorem: for a core S, target
 * g in (0,1/2), integer x not in S,
 *
 *     gap(S+{x}) > g  <=>  not (A)
 *     gap(S+{x}) = g  <=>  (A) and (B)
 *     gap(S+{x}) < g  <=>  (A) and not (B)      [silent killer]
 *
 * (A) kill:  every escape interval (alpha,beta) of S at level g has an
 *            integer m with m-g <= x*alpha and x*beta <= m+g;
 * (B) touch: some level time t in L(S,g) = {g_S = g} has ||x t|| >= g.
 *
 * A silent killer is an escape-killer that touches nowhere: it drops
 * the gap of the extension strictly below g.  At the LRC level
 * g = 1/(m+1), "no silent killers for any (m-1)-core with gap > g" is
 * exactly the induction step LRC(m-1) -> LRC(m) (Theorem R in the
 * silent-killers audit).
 *
 * This program certifies, over ALL (m-1)-subsets S of [1..B] (cores of
 * the m-speed sets), that D(S,g) = X(S,g): every escape-killer touches.
 * Consequence: no m-set T with T \ {max(T)} inside [1..B] has
 * gap(T) < 1/(m+1) (unless a core is itself on/below the level, or a
 * silent killer with x above the radius cap -- both are flagged).
 *
 * Method per core:
 *   probe  midpoints of the crossing grid of the largest speed s_max
 *          (crossings of ||s_max t|| = g, which for g < 1/2 are
 *          generated already sorted).  Any escape of length > 2g
 *          contains a grid point (grid gap (gden-2gnum)/(gden s_max)
 *          < 2g since s_max >= K), so if L_max > 2g some probe has
 *          g_S(t0) > g.  For such a probe the exact local escape
 *          around t0 (nearest crossings of each speed) is computed;
 *          if its length > 2g exactly, then R = floor(2g/L_max) = 0,
 *          D is empty -- core SKIPPED with a certificate.
 *   full   (only near-tight cores: L_max <= 2g) exact event scan:
 *          events = all crossings of all speeds (sorted, deduped);
 *          segment signs at exact midpoints; escapes = maximal
 *          positive runs; level times = events with g_S = g; then D
 *          and X by (A)/(B) for every integer x in [1, min(R, cap)]\S.
 *
 * Modes:
 *   census:  ./lrc_nsk_census <m> <B> <prefix> <nshards> <shard> [cap]
 *            (g = 1/(m+1); cores = (m-1)-subsets of [1..B])
 *   check:   ./lrc_nsk_census check <gnum> <gden> <xmax>
 *            (cores on stdin: "k v1 ... vk"; full exact computation;
 *            prints R, D, X per core -- the Python cross-check port)
 *
 * All arithmetic is exact (int64 rationals); magnitudes stay < 1e9.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;

#define MAXK 12
#define MAXEV 8192
#define MAXCR 8192

static ll gnum, gden;              /* target level g = gnum/gden < 1/2 */
static int K;                      /* core size (m-1) */
static ll S[MAXK];                 /* core speeds, ascending */
static ll cap = 3000;              /* radius cap */

/* ---------- exact integer helpers ---------- */

static ll fdiv_ll(ll n, ll d) {    /* floor(n/d), d > 0 */
    ll q = n / d, r = n % d;
    if (r < 0) q -= 1;
    return q;
}
static ll cdiv_ll(ll n, ll d) {    /* ceil(n/d), d > 0 */
    return -fdiv_ll(-n, d);
}
static ll dist_mod(ll a, ll d) {   /* dist(a mod d, {0,d}), d > 0, a >= 0 */
    ll r = a % d;
    return r < d - r ? r : d - r;
}

/* ||s t|| > g and ||s t|| = g for t = n/d (exact) */
static int saw_gt(ll s, ll n, ll d) {
    return dist_mod(s * n % d, d) * gden > gnum * d;
}
static int saw_eq(ll s, ll n, ll d) {
    return dist_mod(s * n % d, d) * gden == gnum * d;
}
static int gS_gt(ll n, ll d) {     /* g_S(t) > g */
    for (int i = 0; i < K; i++) if (!saw_gt(S[i], n, d)) return 0;
    return 1;
}
static int gS_eq(ll n, ll d) {     /* g_S(t) = g exactly */
    int one = 0;
    for (int i = 0; i < K; i++) {
        if (saw_eq(S[i], n, d)) one = 1;
        else if (!saw_gt(S[i], n, d)) return 0;   /* below g */
    }
    return one;
}

/* ---------- crossings of one speed (generated in sorted order) ----------
 * ||s t|| = g at t = (k*gden - gnum)/(s*gden) and (k*gden+gnum)/(s*gden),
 * k = 0..s, clipped to [0,1].  For 2*gnum < gden the minus/plus forms
 * interleave strictly, so k-order generation is sorted. */
typedef struct { ll num, den; } Rat;

static int rat_cmp(const void *a, const void *b) {
    const Rat *x = (const Rat *)a, *y = (const Rat *)b;
    ll l = x->num * y->den, r = y->num * x->den;
    return l < r ? -1 : (l > r ? 1 : 0);
}

static int crossings_of(ll s, Rat *out) {
    int n = 0;
    ll den = s * gden;
    if (2 * gnum >= gden) {        /* not sorted-generation safe: sort */
        for (ll k = 0; k <= s; k++) {
            ll a = k * gden - gnum, b = k * gden + gnum;
            if (a >= 0 && a <= den && n < MAXCR) { out[n].num = a; out[n].den = den; n++; }
            if (b >= 0 && b <= den && n < MAXCR) { out[n].num = b; out[n].den = den; n++; }
        }
        qsort(out, n, sizeof(Rat), rat_cmp);
        return n;
    }
    for (ll k = 0; k <= s; k++) {
        ll a = k * gden - gnum, b = k * gden + gnum;
        if (a >= 0 && a <= den) { out[n].num = a; out[n].den = den; n++; }
        if (b >= 0 && b <= den) { out[n].num = b; out[n].den = den; n++; }
        if (n >= MAXCR) return 0;  /* guard */
    }
    return n;
}

/* ---------- probe phase: certify L_max > 2g ---------- */

/* nearest crossing of speed s strictly below t0 = n0/d0 (in [0,1]) */
static int nearest_below(ll s, ll n0, ll d0, ll *cnum, ll *cden) {
    ll den = s * gden, D = gden * d0, best = -1;
    ll A = n0 * den + gnum * d0;             /* minus form: k*D < A */
    if (A > 0) {
        ll k = (A - 1) / D;
        if (k <= s) { ll v = k * gden - gnum; if (v >= 0 && v > best) best = v; }
    }
    ll A2 = n0 * den - gnum * d0;            /* plus form: k*D < A2 */
    if (A2 > 0) {
        ll k = (A2 - 1) / D;
        if (k <= s) { ll v = k * gden + gnum; if (v > best) best = v; }
    }
    if (best < 0) return 0;
    *cnum = best; *cden = den;
    return best * d0 < n0 * den;
}
/* nearest crossing of speed s strictly above t0 (in [0,1]) */
static int nearest_above(ll s, ll n0, ll d0, ll *cnum, ll *cden) {
    ll den = s * gden, D = gden * d0, best = -1;
    ll A = n0 * den + gnum * d0;             /* minus form: k*D > A */
    ll k = (A >= 0) ? A / D + 1 : 0;
    if (k <= s) { ll v = k * gden - gnum; if (v >= 0) best = v; }
    ll A2 = n0 * den - gnum * d0;            /* plus form: k*D > A2 */
    ll k2 = (A2 >= 0) ? A2 / D + 1 : 0;
    if (k2 <= s) { ll v = k2 * gden + gnum; if (best < 0 || v < best) best = v; }
    if (best < 0 || best > den) return 0;
    *cnum = best; *cden = den;
    return best * d0 > n0 * den;
}

/* exact local escape around t0 (requires g_S(t0) > g): the escape
 * containing t0, with boundaries at the nearest crossings. */
static int local_escape(ll n0, ll d0, ll *lnum, ll *lden, ll *rnum, ll *rden) {
    ll bn = -1, bd = 1, an = -1, ad = 1;
    for (int i = 0; i < K; i++) {
        ll cn, cd;
        if (!nearest_below(S[i], n0, d0, &cn, &cd)) return 0;
        if (bn < 0 || cn * bd > bn * cd) { bn = cn; bd = cd; }
        if (!nearest_above(S[i], n0, d0, &cn, &cd)) return 0;
        if (an < 0 || cn * ad < an * cd) { an = cn; ad = cd; }
    }
    *lnum = bn; *lden = bd; *rnum = an; *rden = ad;
    return 1;
}

static int kill_one(ll x, Rat a, Rat b);   /* fwd: integer in [x*b-g, x*a+g] */

static int probe_certifies_skip(void) {
    /* Enhanced probe pass: walk the crossing-grid midpoints of s_max;
     * collect the distinct escapes found (each probe inside an escape
     * recovers exactly that escape via nearest crossings); after each
     * new escape, test whether ANY x in [1, R_found] \ S kills ALL
     * found escapes.  D must kill all escapes (found and unfound) and
     * D <= R <= R_found, so an empty found-killer set certifies
     * D = empty.  Returns 1 with that certificate. */
    static Rat cr[MAXCR];
    static Rat found[64][2];
    int nf = 0;
    ll LN = 0, LD = 1;                     /* longest found escape */
    ll smax = S[K - 1];
    int nc = crossings_of(smax, cr);
    if (nc < 2) return 0;
    for (int i = 0; i + 1 < nc; i++) {
        ll d0 = 2 * cr[i].den;
        ll n0 = cr[i].num + cr[i + 1].num;
        if (!gS_gt(n0, d0)) continue;
        ll ln, ld, rn, rd;
        if (!local_escape(n0, d0, &ln, &ld, &rn, &rd)) continue;
        /* dedupe: same escape as one already found? */
        int dup = 0;
        for (int j = 0; j < nf; j++)
            if (found[j][0].num * ld == ln * found[j][0].den &&
                found[j][1].num * rd == rn * found[j][1].den) { dup = 1; break; }
        if (dup) continue;
        if (nf >= 64) return 0;             /* too many: full scan */
        found[nf][0].num = ln; found[nf][0].den = ld;
        found[nf][1].num = rn; found[nf][1].den = rd;
        nf++;
        ll nn = rn * ld - ln * rd, nd = ld * rd;
        if (nn * LD > LN * nd) { LN = nn; LD = nd; }
        /* L_found > 2g  =>  R_found = 0  =>  D empty */
        if (LN * gden > 2 * gnum * LD) return 1;
        /* candidate test: any x <= R_found killing all found escapes? */
        ll Rf = (2 * gnum * LD) / (gden * LN);
        if (Rf > cap) continue;             /* cannot certify cheaply */
        int any = 0;
        for (ll x = 1; x <= Rf && !any; x++) {
            int inS = 0;
            for (int j = 0; j < K; j++) if (S[j] == x) { inS = 1; break; }
            if (inS) continue;
            int ok = 1;
            for (int j = 0; j < nf && ok; j++)
                if (!kill_one(x, found[j][0], found[j][1])) ok = 0;
            if (ok) any = 1;
        }
        if (!any) return 1;                /* D = empty certified */
    }
    return 0;
}

/* ---------- full exact computation ---------- */

typedef struct {
    Rat ev[MAXEV];      /* sorted unique events (crossings) */
    int nev;
    Rat esc[2 * MAXEV]; /* escapes as (a,b) pairs */
    int nesc;
    Rat lvl[MAXEV];     /* level times */
    int nlvl;
} CoreData;

static int full_scan(CoreData *cd) {
    static Rat all[MAXEV];
    int n = 0;
    for (int i = 0; i < K; i++) {
        Rat tmp[MAXCR];
        int c = crossings_of(S[i], tmp);
        if (c == 0) return -1;
        if (n + c > MAXEV) return -1;
        memcpy(all + n, tmp, c * sizeof(Rat));
        n += c;
    }
    qsort(all, n, sizeof(Rat), rat_cmp);
    int nev = 0;
    for (int i = 0; i < n; i++) {
        if (nev && all[i].num * cd->ev[nev-1].den == cd->ev[nev-1].num * all[i].den)
            continue;
        cd->ev[nev++] = all[i];
    }
    cd->nev = nev;

    static char pos[MAXEV];
    for (int i = 0; i + 1 < nev; i++) {
        ll dm = 2 * cd->ev[i].den * cd->ev[i+1].den;
        ll nm = cd->ev[i].num * cd->ev[i+1].den
              + cd->ev[i+1].num * cd->ev[i].den;
        pos[i] = (char)gS_gt(nm, dm);
    }
    cd->nlvl = 0;
    for (int i = 0; i < nev; i++)
        if (gS_eq(cd->ev[i].num, cd->ev[i].den))
            cd->lvl[cd->nlvl++] = cd->ev[i];

    cd->nesc = 0;
    int i = 0;
    while (i + 1 < nev) {
        if (!pos[i]) { i++; continue; }
        int j = i;
        while (j + 1 < nev && pos[j + 1]
               && gS_gt(cd->ev[j+1].num, cd->ev[j+1].den))
            j++;
        cd->esc[2 * cd->nesc]     = cd->ev[i];
        cd->esc[2 * cd->nesc + 1] = cd->ev[j + 1];
        cd->nesc++;
        i = j + 1;
    }
    return cd->nesc;
}

/* ---------- D and X (kill_one declared before probe) ---------- */

static int kill_one(ll x, Rat a, Rat b) {   /* integer in [x*b-g, x*a+g] */
    ll nb = x * b.num * gden - gnum * b.den, db = b.den * gden;
    ll na = x * a.num * gden + gnum * a.den, da = a.den * gden;
    return cdiv_ll(nb, db) <= fdiv_ll(na, da);
}
static int kill_all(ll x, CoreData *cd) {
    for (int j = 0; j < cd->nesc; j++)
        if (!kill_one(x, cd->esc[2*j], cd->esc[2*j+1])) return 0;
    return 1;
}
static int touch(ll x, CoreData *cd) {
    for (int j = 0; j < cd->nlvl; j++) {
        Rat t = cd->lvl[j];
        if (dist_mod(x * t.num % t.den, t.den) * gden >= gnum * t.den)
            return 1;
    }
    return 0;
}

/* ---------- census ---------- */

static FILE *f_cores, *f_flags;
static ll n_total, n_skip, n_full, n_nonemptyD, n_silent, n_on, n_sub,
          n_rflag, max_R, sumD, sumX, max_silent_x;

static void print_S(FILE *f) {
    fprintf(f, "S=");
    for (int i = 0; i < K; i++) fprintf(f, "%lld,", S[i]);
}

static void process_core(void) {
    n_total++;
    if (probe_certifies_skip()) { n_skip++; return; }

    CoreData cd;
    memset(&cd, 0, sizeof cd);
    int nesc = full_scan(&cd);
    if (nesc < 0) {                   /* overflow guard: treat as full */
        fprintf(f_flags, "OVF  ");
        print_S(f_flags);
        fprintf(f_flags, "\n");
        return;
    }
    if (nesc == 0) {                  /* gap(S) <= g */
        if (cd.nlvl > 0) n_on++; else n_sub++;
        fprintf(f_flags, cd.nlvl ? "ON   " : "SUB  ");
        print_S(f_flags);
        fprintf(f_flags, " (gap %s %lld/%lld at core level %d)\n",
                cd.nlvl ? "=" : "<", gnum, gden, K);
        return;
    }
    n_full++;

    ll LN = 0, LD = 1;
    for (int j = 0; j < nesc; j++) {
        Rat a = cd.esc[2*j], b = cd.esc[2*j+1];
        ll nn = b.num * a.den - a.num * b.den, nd = a.den * b.den;
        if (nn * LD > LN * nd) { LN = nn; LD = nd; }
    }
    ll R = (2 * gnum * LD) / (gden * LN);         /* floor(2g / L_max) */
    if (R > max_R) max_R = R;
    ll xmax = R < cap ? R : cap;
    if (R > cap) {
        n_rflag++;
        fprintf(f_flags, "RFLAG ");
        print_S(f_flags);
        fprintf(f_flags, "  R=%lld > cap=%lld\n", R, cap);
    }

    ll Dbuf[512], Xbuf[512];
    int nD = 0, nX = 0;
    for (ll x = 1; x <= xmax; x++) {
        int inS = 0;
        for (int i = 0; i < K; i++) if (S[i] == x) { inS = 1; break; }
        if (inS || !kill_all(x, &cd)) continue;
        if (nD < 512) Dbuf[nD] = x;
        nD++;
        if (touch(x, &cd)) {
            if (nX < 512) Xbuf[nX] = x;
            nX++;
        } else {
            n_silent++;
            if (x > max_silent_x) max_silent_x = x;
            fprintf(f_flags, "SILENT ");
            print_S(f_flags);
            fprintf(f_flags, "  x=%lld R=%lld\n", x, R);
        }
    }
    if (nD > 0) {
        n_nonemptyD++;
        sumD += nD; sumX += nX;
        print_S(f_cores);
        fprintf(f_cores, " | R=%lld | D=%lld", R, Dbuf[0]);
        for (int i = 1; i < nD && i < 512; i++) fprintf(f_cores, ",%lld", Dbuf[i]);
        fprintf(f_cores, " | X=%lld", nX ? Xbuf[0] : -1LL);
        for (int i = 1; i < nX && i < 512; i++) fprintf(f_cores, ",%lld", Xbuf[i]);
        fprintf(f_cores, "\n");
        if (nD > 512 || nX > 512)
            fprintf(f_flags, "LONG  "), print_S(f_flags),
            fprintf(f_flags, " D=%d X=%d truncated\n", nD, nX);
    }
}

static ll idx_counter, shard, nshards;

static void recurse(int depth, ll start, ll B) {
    if (depth == K) {
        if (idx_counter % nshards == shard) process_core();
        idx_counter++;
        return;
    }
    for (ll v = start; v <= B; v++) {
        S[depth] = v;
        recurse(depth + 1, v + 1, B);
    }
}

/* ---------- check mode ---------- */

static int check_mode(ll gn, ll gd, ll xmax) {
    gnum = gn; gden = gd;
    char line[8192];
    while (fgets(line, sizeof line, stdin)) {
        int k;
        if (sscanf(line, "%d", &k) != 1 || k < 1 || k > MAXK) continue;
        char *p = line;
        while (*p == ' ') p++;
        while (*p && *p != ' ') p++;
        ll v[MAXK];
        int ok = 1;
        for (int i = 0; i < k; i++) {
            if (sscanf(p, "%lld", &v[i]) != 1) { ok = 0; break; }
            while (*p == ' ') p++;
            while (*p && *p != ' ') p++;
        }
        if (!ok) continue;
        K = k;
        for (int i = 0; i < K; i++) S[i] = v[i];
        CoreData cd;
        memset(&cd, 0, sizeof cd);
        int nesc = full_scan(&cd);
        if (nesc < 0) { printf("OVERFLOW\n"); continue; }
        print_S(stdout);
        if (nesc == 0) {
            printf(" | %s\n", cd.nlvl ? "ON-level (gap=g)" : "BELOW (gap<g)");
            continue;
        }
        ll LN = 0, LD = 1;
        for (int j = 0; j < nesc; j++) {
            Rat a = cd.esc[2*j], b = cd.esc[2*j+1];
            ll nn = b.num * a.den - a.num * b.den, nd = a.den * b.den;
            if (nn * LD > LN * nd) { LN = nn; LD = nd; }
        }
        ll R = (2 * gnum * LD) / (gden * LN);
        ll xm = R < xmax ? R : xmax;
        printf(" | R=%lld | D=", R);
        int first = 1;
        for (ll x = 1; x <= xm; x++) {
            int inS = 0;
            for (int i = 0; i < K; i++) if (S[i] == x) { inS = 1; break; }
            if (inS || !kill_all(x, &cd)) continue;
            printf(first ? "%lld" : ",%lld", x);
            first = 0;
        }
        printf(" | X=");
        first = 1;
        for (ll x = 1; x <= xm; x++) {
            int inS = 0;
            for (int i = 0; i < K; i++) if (S[i] == x) { inS = 1; break; }
            if (inS || !kill_all(x, &cd) || !touch(x, &cd)) continue;
            printf(first ? "%lld" : ",%lld", x);
            first = 0;
        }
        printf("\n");
    }
    return 0;
}

/* ---------- main ---------- */

int main(int argc, char **argv) {
    if (argc >= 2 && !strcmp(argv[1], "check")) {
        if (argc < 5) { fprintf(stderr, "check gnum gden xmax\n"); return 2; }
        return check_mode(atoll(argv[2]), atoll(argv[3]), atoll(argv[4]));
    }
    if (argc < 6) {
        fprintf(stderr,
            "usage: %s <m> <B> <prefix> <nshards> <shard> [cap]\n"
            "       %s check <gnum> <gden> <xmax>\n", argv[0], argv[0]);
        return 2;
    }
    ll m = atoll(argv[1]);
    ll B = atoll(argv[2]);
    const char *prefix = argv[3];
    nshards = atoll(argv[4]);
    shard = atoll(argv[5]);
    if (argc > 6) cap = atoll(argv[6]);
    if (m < 2 || m > 11 || B < 2 || nshards < 1 || shard < 0
        || shard >= nshards) {
        fprintf(stderr, "bad args\n");
        return 2;
    }
    gnum = 1; gden = m + 1;
    K = (int)m - 1;

    char path[1024];
    snprintf(path, sizeof path, "%s_cores.txt", prefix);
    f_cores = fopen(path, "w");
    snprintf(path, sizeof path, "%s_flags.txt", prefix);
    f_flags = fopen(path, "w");
    if (!f_cores || !f_flags) { perror("open"); return 1; }

    recurse(0, 1, B);

    fclose(f_cores); fclose(f_flags);

    snprintf(path, sizeof path, "%s_summary.txt", prefix);
    FILE *fs = fopen(path, "w");
    FILE *out[2] = { stdout, fs };
    for (int o = 0; o < 2; o++) {
        fprintf(out[o],
            "m=%lld B=%lld g=1/%lld cores=%lld shard=%lld/%lld cap=%lld\n",
            m, B, gden, n_total, shard, nshards, cap);
        fprintf(out[o],
            "skip_certified=%lld full=%lld nonemptyD=%lld sumD=%lld sumX=%lld\n",
            n_skip, n_full, n_nonemptyD, sumD, sumX);
        fprintf(out[o],
            "SILENT=%lld max_silent_x=%lld ON=%lld SUB=%lld RFLAG=%lld "
            "max_R=%lld\n",
            n_silent, max_silent_x, n_on, n_sub, n_rflag, max_R);
    }
    fclose(fs);
    return 0;
}
