/* E-004: reimplementa a regra impressa por Volkov (2006), secao 5.1, p.19.
 *
 * O MAPA E' O DA EQUACAO (1.1) DELE, p.2, e nao o acelerado:
 *
 *     M(x) = 5x + 1   se x e' impar
 *     M(x) = x / 2    se x e' par
 *
 * Uma primeira versao deste arquivo usava (5x+1)/2, que e' o T_5 de
 * Kontorovich-Lagarias, um mapa diferente. A orbita do Volkov passa pelo
 * pico 5x+1; aquela versao registrava metade dele, entao o teste do
 * limiar nunca via o pico, e toda semente cujo pico caisse em (T, 2T]
 * era classificada errado. Achado pela rodada 3 de critica.
 *
 * Regra dele: um k <= m e' "good" se a orbita entra num dos tres ciclos
 * positivos conhecidos, e "bad" se excede o limiar T que faz as vezes de
 * infinito. Os ciclos sao calculados, nao copiados: o proprio texto dele
 * imprime o primeiro como 1 -> 6 -> 3 -> 1, mas M(3) = 16, entao o ciclo
 * impresso esta incompleto. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef unsigned long long u64;

static u64 CYC[64]; static int NCYC = 0;

static void build_cycles(void){
    const u64 seeds[3] = {1, 13, 17};
    for (int s = 0; s < 3; s++){
        u64 path[512]; int len = 0; u64 n = seeds[s];
        for (;;){
            int seen = -1;
            for (int i = 0; i < len; i++) if (path[i] == n){ seen = i; break; }
            if (seen >= 0){ for (int i = seen; i < len; i++) CYC[NCYC++] = path[i]; break; }
            path[len++] = n;
            n = (n % 2) ? 5*n + 1 : n / 2;
        }
    }
}
static int in_cycle(u64 n){
    if (n > 416) return 0;
    for (int i = 0; i < NCYC; i++) if (CYC[i] == n) return 1;
    return 0;
}
int main(int argc, char **argv){
    if (argc < 3){ fprintf(stderr, "uso: %s <m> <T>\n", argv[0]); return 2; }
    build_cycles();
    u64 m = strtoull(argv[1], 0, 10), T = strtoull(argv[2], 0, 10);
    u64 good = 0, bad = 0, other = 0;
    for (u64 k = 1; k <= m; k++){
        u64 n = k; int v = 0;
        for (int s = 0; s < 400000; s++){
            if (n > T){ v = -1; break; }
            if (in_cycle(n)){ v = 1; break; }
            n = (n % 2) ? 5*n + 1 : n / 2;
        }
        if (v == 1) good++; else if (v == -1) bad++; else other++;
    }
    printf("m=%llu T=%.3e  good=%llu bad=%llu other=%llu  beta=%.5f\n",
           m, (double)T, good, bad, other,
           good ? log((double)good)/log((double)m) : 0.0);
    return 0;
}
