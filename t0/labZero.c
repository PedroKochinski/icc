#include <float.h>
#include <math.h>
#include <stdio.h>

#include "ZeroFuncao.h"
#include "utils.h"



int main() {
    int it = 0, tipoCalculo = 0;
    real_t a, b;
    Polinomio pol;
    real_t raiz = 0.0;
    real_t erro = 0.0;
    double tempo = 0.0;

    scanf("%d", &pol.grau);  // le o grau do polinomio
    pol.p = (real_t*)malloc((pol.grau + 1) * sizeof(real_t)); // aloca o vetor de coeficientes
    for (int i = pol.grau; i >= 0; --i) {
        scanf("%lf", &pol.p[i]);  // le os coeficientes do polinomio
    }
    scanf("%lf %lf", &a, &b);  // intervalo onde está uma das raizes.

    
    printf("RAPIDO\n\n");
    for (int i = 1; i < 4; i++) {
        it = 0;
        raiz = 0.0;
        tempo = timestamp();
        erro = bisseccao(pol, a, b, i, &it, &raiz, tipoCalculo);
        tempo = timestamp() - tempo;
        printf("bissec <raiz: %.15lf> <criterio(%d): %lf> <it: %d> <tempo: %f>\n", raiz, i, erro, it, tempo);
    }
    for (int i = 1; i < 4; i++) {
        it = 0;
        raiz = 0.0;
        tempo = timestamp();
        erro = newtonRaphson(pol, a, i, &it, &raiz, tipoCalculo);
        tempo = timestamp() - tempo;
        printf("newton <raiz: %.15lf> <criterio(%d): %lf> <it: %d> <tempo: %f>\n", raiz, i, erro, it, tempo);
    }

    tipoCalculo = 1;  // tipoCalculo = 1 -> lento
    printf("\n");
    printf("LENTO\n\n");
    for (int i = 1; i < 4; i++) {
        it = 0;
        raiz = 0.0;
        tempo = timestamp();
        erro = bisseccao(pol, a, b, i, &it, &raiz, tipoCalculo);
        tempo = timestamp() - tempo;
        printf("bissec <raiz: %.15lf> <criterio(%d): %lf> <it: %d> <tempo: %f>\n", raiz, i, erro, it, tempo);
    }

    for (int i = 1; i < 4; i++) {
        it = 0;
        raiz = 0.0;
        tempo = timestamp();
        erro = newtonRaphson(pol, a, i, &it, &raiz, tipoCalculo);
        tempo = timestamp() - tempo;
        printf("newton <raiz: %.15lf> <criterio(%d): %lf> <it: %d> <tempo: %f>\n", raiz, i, erro, it, tempo);
    }
    return 0;
}
