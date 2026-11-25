/*
 * tmr.c - Redundancia Modular Tripla (TMR) em software:
 *         custo de desempenho e cobertura de falhas.
 *
 * Cargas:   mm  = multiplicacao de matrizes NxN (uint32, limitada por CPU)
 *           vs  = out[i] = a*x[i] + b (uint32, limitada por memoria)
 * Esquemas: base (1 replica), dmr (2 replicas + comparacao),
 *           tmr (3 replicas sequenciais + votador bit a bit),
 *           tmrp (3 replicas em threads + votador bit a bit)
 *
 * Uso:
 *   ./tmr time   <mm|vs> <tamanho> <reps>          -> CSV de tempos
 *   ./tmr voter  <n_palavras> <reps>               -> CSV votador palavra x bit a bit
 *   ./tmr inject <mm|vs> <tamanho> <n_inj> <seed>  -> CSV de resultados de injecao
 *
 * Compilar: gcc -O2 -march=native -pthread -o tmr tmr.c
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

typedef uint32_t u32;

/* ---------------- utilidades ---------------- */
static double now(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}
static uint64_t rng_s = 88172645463325252ULL;
static inline uint64_t rng(void) {            /* xorshift64 */
    rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s;
}
static void *xmalloc(size_t n) {
    void *p = aligned_alloc(64, (n + 63) / 64 * 64);
    if (!p) { perror("malloc"); exit(1); }
    return p;
}

/* ---------------- cargas ---------------- */
typedef struct {
    int kind;          /* 0 = mm, 1 = vs */
    size_t n;          /* N da matriz ou tamanho do vetor */
    size_t in_len;     /* palavras de entrada */
    size_t out_len;    /* palavras de saida  */
} work_t;

static void compute(const work_t *w, const u32 *in, u32 *out) {
    if (w->kind == 0) {
        size_t N = w->n; const u32 *A = in, *B = in + N * N;
        memset(out, 0, N * N * sizeof(u32));
        for (size_t i = 0; i < N; i++)
            for (size_t k = 0; k < N; k++) {
                u32 a = A[i * N + k];
                const u32 *Bk = B + k * N; u32 *Ci = out + i * N;
                for (size_t j = 0; j < N; j++) Ci[j] += a * Bk[j];
            }
    } else {
        const u32 a = 2654435761u, b = 12345u;
        for (size_t i = 0; i < w->n; i++) out[i] = a * in[i] + b;
    }
}

static work_t make_work(const char *k, size_t n) {
    work_t w; w.n = n;
    if (!strcmp(k, "mm")) { w.kind = 0; w.in_len = 2 * n * n; w.out_len = n * n; }
    else                  { w.kind = 1; w.in_len = n;         w.out_len = n; }
    return w;
}

/* ---------------- votadores ---------------- */
/* bit a bit: V = AB + AC + BC ; D = (A^B) | (A^C)  (sem desvios) */
static int vote_bitwise(const u32 *a, const u32 *b, const u32 *c, u32 *v, size_t n) {
    u32 d = 0;
    for (size_t i = 0; i < n; i++) {
        u32 x = a[i], y = b[i], z = c[i];
        v[i] = (x & y) | (x & z) | (y & z);
        d |= (x ^ y) | (x ^ z);
    }
    return d != 0;                      /* 1 = houve divergencia */
}
/* por palavra: escolhe o valor que aparece >= 2 vezes; se os tres
 * forem diferentes, nao ha maioria (erro detectado nao corrigivel). */
static int vote_word(const u32 *a, const u32 *b, const u32 *c, u32 *v, size_t n, int *nomaj) {
    int d = 0; *nomaj = 0;
    for (size_t i = 0; i < n; i++) {
        u32 x = a[i], y = b[i], z = c[i];
        if (x == y || x == z) { v[i] = x; if (x != y || x != z) d = 1; }
        else if (y == z)      { v[i] = y; d = 1; }
        else                  { v[i] = x; d = 1; *nomaj = 1; }
    }
    return d;
}
static int compare(const u32 *a, const u32 *b, size_t n) {
    u32 d = 0;
    for (size_t i = 0; i < n; i++) d |= a[i] ^ b[i];
    return d != 0;
}

/* ---------------- replicas ---------------- */
typedef struct {
    const work_t *w; const u32 *master; u32 *in; u32 *out;
} rep_t;
static void run_replica(rep_t *r) {
    memcpy(r->in, r->master, r->w->in_len * sizeof(u32));   /* copia privada */
    compute(r->w, r->in, r->out);
}
static void *thr(void *p) { run_replica((rep_t *)p); return NULL; }

/* ---------------- modo time ---------------- */
static void mode_time(const char *k, size_t n, int reps) {
    work_t w = make_work(k, n);
    u32 *master = xmalloc(w.in_len * 4);
    for (size_t i = 0; i < w.in_len; i++) master[i] = (u32)rng();
    u32 *in[3], *out[3], *v = xmalloc(w.out_len * 4);
    for (int r = 0; r < 3; r++) { in[r] = xmalloc(w.in_len * 4); out[r] = xmalloc(w.out_len * 4);
        memset(in[r], 0, w.in_len * 4); memset(out[r], 0, w.out_len * 4); }
    memset(v, 0, w.out_len * 4);
    rep_t R[3]; for (int r = 0; r < 3; r++) R[r] = (rep_t){&w, master, in[r], out[r]};

    const char *sch[] = {"base", "dmr", "tmr", "tmrp"};
    printf("carga,tamanho,esquema,rep,t_total,t_replicas,t_voto\n");
    for (int warm = 0; warm < 2; warm++) run_replica(&R[0]);       /* aquecimento */
    for (int rep = 0; rep < reps; rep++)
        for (int s = 0; s < 4; s++) {                               /* intercalado */
            double t0 = now(), t1, t2;
            if (s == 0) { run_replica(&R[0]); t1 = now(); }
            else if (s == 1) { run_replica(&R[0]); run_replica(&R[1]); t1 = now();
                               v[0] = (u32)compare(out[0], out[1], w.out_len); }
            else if (s == 2) { for (int r = 0; r < 3; r++) run_replica(&R[r]); t1 = now();
                               vote_bitwise(out[0], out[1], out[2], v, w.out_len); }
            else { pthread_t th[3];
                   for (int r = 0; r < 3; r++) pthread_create(&th[r], NULL, thr, &R[r]);
                   for (int r = 0; r < 3; r++) pthread_join(th[r], NULL);
                   t1 = now(); vote_bitwise(out[0], out[1], out[2], v, w.out_len); }
            t2 = now();
            printf("%s,%zu,%s,%d,%.6f,%.6f,%.6f\n", k, n, sch[s], rep, t2 - t0, t1 - t0, t2 - t1);
        }
    volatile u32 sink = v[w.out_len / 2]; (void)sink;
}

/* ---------------- modo voter ---------------- */
static void mode_voter(size_t n, int reps) {
    u32 *a = xmalloc(n * 4), *b = xmalloc(n * 4), *c = xmalloc(n * 4), *v = xmalloc(n * 4);
    printf("cenario,votador,rep,ns_por_palavra\n");
    /* cenario "limpo": replicas iguais; "ruidoso": 50% das palavras com 1 replica divergente
     * em posicao aleatoria (estressa o preditor de desvios). */
    for (int scen = 0; scen < 2; scen++) {
        for (size_t i = 0; i < n; i++) { a[i] = b[i] = c[i] = (u32)rng(); }
        if (scen == 1) for (size_t i = 0; i < n; i++) if (rng() & 1) {
            int r = rng() % 3; u32 f = 1u << (rng() % 32);
            if (r == 0) a[i] ^= f; else if (r == 1) b[i] ^= f; else c[i] ^= f; }
        for (int rep = 0; rep < reps; rep++) for (int vt = 0; vt < 2; vt++) {
            int nm; double t0 = now();
            if (vt == 0) vote_word(a, b, c, v, n, &nm); else vote_bitwise(a, b, c, v, n);
            double t = now() - t0;
            printf("%s,%s,%d,%.4f\n", scen ? "ruidoso" : "limpo", vt ? "bit_a_bit" : "palavra", rep, t * 1e9 / n);
        }
    }
    volatile u32 sink = v[n / 2]; (void)sink;
}

/* ---------------- modo inject ---------------- */
/* Cenarios (modelo de falha: inversao de bit unico, SEU):
 *  in1   : 1 bit na entrada privada de UMA replica (antes do calculo)
 *  out1  : 1 bit na saida de UMA replica (antes do voto)
 *  out2d : 2 replicas, palavras de saida diferentes
 *  out2w : 2 replicas, mesma palavra, bits diferentes
 *  out2s : 2 replicas, mesma palavra, mesmo bit (erro identico)
 *  comum : 1 bit na entrada mestre (modo comum, antes da replicacao)
 * Classes: OK (saida correta, sem alarme), CORR (corrigida, com alarme D),
 *          DUE (detectado, nao corrigido -> sistema para/reexecuta),
 *          SDC (saida errada sem alarme ou com voto errado). */
enum { OK, CORR, DUE, SDC };
static const char *cls_name[] = {"OK", "CORR", "DUE", "SDC"};

static void mode_inject(const char *k, size_t n, int ninj, uint64_t seed) {
    rng_s = seed;
    work_t w = make_work(k, n);
    u32 *master = xmalloc(w.in_len * 4), *gold = xmalloc(w.out_len * 4), *v = xmalloc(w.out_len * 4);
    for (size_t i = 0; i < w.in_len; i++) master[i] = (u32)rng();
    compute(&w, master, gold);
    u32 *mc = xmalloc(w.in_len * 4);
    u32 *in[3], *out[3];
    for (int r = 0; r < 3; r++) { in[r] = xmalloc(w.in_len * 4); out[r] = xmalloc(w.out_len * 4); }

    const char *scen[] = {"in1", "out1", "out2d", "out2w", "out2s", "comum"};
    const char *sch[] = {"base", "dmr", "tmr_palavra", "tmr_bit"};
    int cnt[6][4][4]; memset(cnt, 0, sizeof cnt);

    for (int s = 0; s < 6; s++)
        for (int it = 0; it < ninj; it++) {
            /* sorteio da falha (o mesmo para todos os esquemas) */
            size_t ei = rng() % w.in_len, eo = rng() % w.out_len, eo2 = rng() % w.out_len;
            while (eo2 == eo) eo2 = rng() % w.out_len;
            int b1 = rng() % 32, b2 = rng() % 32; while (b2 == b1) b2 = rng() % 32;
            int ra = rng() % 3, rb = (ra + 1 + rng() % 2) % 3;   /* replicas atingidas */

            memcpy(mc, master, w.in_len * 4);
            if (s == 5) mc[ei] ^= 1u << b1;                        /* modo comum */
            for (int r = 0; r < 3; r++) {
                memcpy(in[r], mc, w.in_len * 4);
                if (s == 0 && r == ra) in[r][ei] ^= 1u << b1;
                compute(&w, in[r], out[r]);
            }
            if (s == 1) out[ra][eo] ^= 1u << b1;
            if (s == 2) { out[ra][eo] ^= 1u << b1; out[rb][eo2] ^= 1u << b2; }
            if (s == 3) { out[ra][eo] ^= 1u << b1; out[rb][eo]  ^= 1u << b2; }
            if (s == 4) { out[ra][eo] ^= 1u << b1; out[rb][eo]  ^= 1u << b1; }

            /* base: usa a replica ra (a atingida, para falhas de 1 replica);
             * dmr: replicas {ra, rb}; tmr: as tres. */
            for (int c = 0; c < 4; c++) {
                int cl;
                if (c == 0) {
                    cl = memcmp(out[ra], gold, w.out_len * 4) ? SDC : OK;
                } else if (c == 1) {
                    if (compare(out[ra], out[rb], w.out_len)) cl = DUE;
                    else cl = memcmp(out[ra], gold, w.out_len * 4) ? SDC : OK;
                } else {
                    int d, nm = 0;
                    if (c == 2) d = vote_word(out[0], out[1], out[2], v, w.out_len, &nm);
                    else        d = vote_bitwise(out[0], out[1], out[2], v, w.out_len);
                    int wrong = memcmp(v, gold, w.out_len * 4) != 0;
                    if (nm) cl = DUE;
                    else if (wrong) cl = SDC;
                    else cl = d ? CORR : OK;
                }
                cnt[s][c][cl]++;
            }
        }
    printf("carga,tamanho,cenario,esquema,OK,CORR,DUE,SDC\n");
    for (int s = 0; s < 6; s++) for (int c = 0; c < 4; c++)
        printf("%s,%zu,%s,%s,%d,%d,%d,%d\n", k, n, scen[s], sch[c],
               cnt[s][c][OK], cnt[s][c][CORR], cnt[s][c][DUE], cnt[s][c][SDC]);
    (void)cls_name;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "uso: ver cabecalho\n"); return 1; }
    if (!strcmp(argv[1], "time") && argc == 5) mode_time(argv[2], strtoull(argv[3], 0, 10), atoi(argv[4]));
    else if (!strcmp(argv[1], "voter") && argc == 4) mode_voter(strtoull(argv[2], 0, 10), atoi(argv[3]));
    else if (!strcmp(argv[1], "inject") && argc == 6)
        mode_inject(argv[2], strtoull(argv[3], 0, 10), atoi(argv[4]), strtoull(argv[5], 0, 10));
    else { fprintf(stderr, "argumentos invalidos\n"); return 1; }
    return 0;
}
