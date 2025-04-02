#include "SistemaLinear.h"

#include <stdio.h>
#include <stdlib.h>

void le_s_linear(s_linear *sistema) {
    for (int i = 0; i < sistema->tamanho; ++i) {
        for (int j = 0; j < sistema->tamanho; j++) {
            scanf("%lf", sistema->A[i][j]);
        }
    }
}

s_linear *cria_s_linear(int tamanho) {
    printf("CRIANDO SLINEAR DE TAMANHO %d\n", tamanho);
    s_linear *s = malloc(sizeof(s_linear*));
    s->A = (real_t **)malloc((tamanho) * sizeof(real_t *));
    // aloca um vetor com todos os elementos da matriz
    s->A[0] = malloc(tamanho * tamanho * sizeof(real_t));
    // ajusta os demais ponteiros de linhas (i > 0)
    for (int i = 1; i < tamanho; i++) s->A[i] = s->A[0] + i * tamanho;
    return s;
}
