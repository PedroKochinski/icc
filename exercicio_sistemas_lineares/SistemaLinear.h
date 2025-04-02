#ifndef _SISTEMALINEAR_H_
#define _SISTEMALINEAR_H_

#include "utils.h"

typedef struct {
    // AX = B
    int tamanho;
    real_t **A;
    real_t *X;
    real_t *B;
} s_linear;

void le_s_linear(s_linear *sistema);

s_linear *cria_s_linear(int tamanho);

#endif