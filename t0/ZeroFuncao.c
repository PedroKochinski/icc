#include "ZeroFuncao.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "utils.h"


// Converte double para int64 preservando os bits
int64_t double_to_int64(double d) {
    int64_t i;
    memcpy(&i, &d, sizeof(double));
    // Ajuste para que números negativos fiquem ordenados corretamente
    return (i < 0) ? INT64_MIN - i : i;
}

int64_t calcula_diferenca_ulp(double a, double b) {
    // Se forem exatamente iguais, a diferença é zero
    if (a == b) return 0;

    int64_t ia = double_to_int64(a);
    int64_t ib = double_to_int64(b);
    int64_t diff = llabs(ia - ib);
    return diff;
}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson(Polinomio p, real_t x0, int criterioParada, int *it,
                     real_t *raiz) {
    int parada = 0;

    real_t derivada = 0.1, ultimaDerivada, ultimaRaiz, novoX = x0, erro = 0.0, px;
    
    while (!parada) {
        (*it)++;
        ultimaRaiz = *raiz;
        ultimaDerivada = derivada;
        px = 0;
        calcPolinomio_rapido(p, *raiz, &px, &derivada);
        // printf("raiz candidata: %lf\n", *raiz);
        
        
        if(derivada <= DBL_EPSILON) {
            derivada = ultimaDerivada;
            // printf("Derivada muito pequena, usando a derivada anterior: %lf\n", derivada);
        }
        
        x0 = *raiz;
        *raiz = *raiz - (px / derivada);  // calcula a raiz candidata
        erro = fabs(*raiz - x0) / *raiz * 100;

        if (*it >= MAXIT) {
            parada = 1;
            // printf("Número máximo de iterações atingido\n");
            return erro;  // iterações excedidas, retorna o erro
        }
        if ((criterioParada == 1) && (*it > 2) && (fabs(ultimaRaiz - (*raiz)) < 0.000001)) {  // criterio 1 - erro absoluto entre xk e xk-1
            parada = 1;
            // printf("Erro absoluto entre xk (%lf) e xk-1 (%lf): %lf\n", ultimaRaiz, *raiz, fabs(ultimaRaiz - (*raiz)));
            return erro;
        }
        if (criterioParada == 2 && fabs(*raiz) <= DBL_EPSILON) {  // criterio 2 - se o valor da de p(x) é quae zero, i.e |p(x)| <= epsilon
            parada = 1;
            // printf("Valor de p(x) é zero: %lf\n", fabs(calcRaizCandidata));
            return erro;
        }
        if (criterioParada == 3 && (*it > 2) && calcula_diferenca_ulp(ultimaRaiz, *raiz) <= 2) {  // criterio 3 - ULP's entre xk e xk-1 <= 2
            parada = 1;
            return erro;
        }

        
    }
}

/*
Polinomio p - polinomio de entrada
real_t a - primeiro valor do intervalo
real_t b - segundo valor do intervalo
int criterioParada - numero maximo de iterações antes de finalizar
int *it - quantidade de iteraçções
real_t raiz - valor da raiz que foi encontrado

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
*/
real_t bisseccao(Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz) {
    int parada = 0;
    real_t erro = 0;
    real_t calcRaizCandidata, calcA, calcB;
    real_t ultimaRaiz = *raiz;
    while (!parada) {
        ultimaRaiz = *raiz;
        *raiz = (a + b) / 2;  // calcula a raiz candidata, o ponto medio do intervalo
        erro = (fabs(ultimaRaiz - *raiz) / ultimaRaiz) * 100;
        // printf("testando a raiz candidata: %lf\n", *raiz);
        calcPolinomio_rapido(p, a, &calcA, NULL);
        calcPolinomio_rapido(p, b, &calcB, NULL);
        calcPolinomio_rapido(p, *raiz, &calcRaizCandidata, NULL);

        if (calcA * calcRaizCandidata < 0) {
            // printf("b recebe raizCandidata, pois a raiz está no intervalo [a, raiz]\n");
            b = *raiz;  // B recebe raizCandidata, pois a raiz está no intervalo [a, raiz]
        } else if (calcA * calcRaizCandidata > 0) {
            // printf("a recebe raizCandidata, pois a raiz está no intervalo [raiz, b]\n");
            a = *raiz;  // A recebe raizCandidata, pois a raiz está no intervalo [raiz, b]
        } else {
            // printf("raiz encontrada\n");
            parada = 1;
            return erro;  // achou a raiz, retorna o erro
        }

        (*it)++;
        // printf("Iteração %d\n", *it);

        if (*it >= MAXIT) {
            parada = 1;
            // printf("Número máximo de iterações atingido\n");
            return erro;  // iterações excedidas, retorna o erro
        }
        if ((criterioParada == 1) && (*it > 2) && (fabs(ultimaRaiz - (*raiz)) < 0.000001)) {  // criterio 1 - erro absoluto entre xk e xk-1
            parada = 1;
            // printf("Erro absoluto entre xk (%lf) e xk-1 (%lf): %lf\n", ultimaRaiz, *raiz, fabs(ultimaRaiz - (*raiz)));
            return erro;
        }
        if (criterioParada == 2 && fabs(calcRaizCandidata) <= DBL_EPSILON) {  // criterio 2 - se o valor da de p(x) é quae zero, i.e |p(x)| <= epsilon
            parada = 1;
            // printf("Valor de p(x) é zero: %lf\n", fabs(calcRaizCandidata));
            return erro;
        }
        if (criterioParada == 3 && (*it > 2) && calcula_diferenca_ulp(ultimaRaiz, *raiz) <= 2) {  // criterio 3 - ULP's entre xk e xk-1 <= 2
            parada = 1;
            return erro;
        }
    }

    return erro;
}

/*
Polinomio p - polinomio de entrada
real_t x - x de entrada
real_t *px - valor de P(X) calculado
real_t *dpx - valor da derivada de P(x) calculado, i.e, f'(x)
*/
void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx) {
    real_t calculo_px = 0;
    real_t calculo_dpx = 0;
    for (int i = p.grau; i > 0; i--) {
        calculo_px = calculo_px * x + p.p[i];
        calculo_dpx = calculo_dpx * x + calculo_px;
    }
    calculo_px = calculo_px * x + p.p[0];
    *px = calculo_px;
    if (dpx != NULL) *dpx = calculo_dpx;
}

/*
Polinomio p - polinomio de entrada
real_t x - x de entrada
real_t *px - valor de P(X) calculado
real_t *dpx - valor da derivada de P(x) calculado, i.e, f'(x)
*/
void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx) {
    for (int i = p.grau; i >= 0; i--) {
        *px += p.p[i] * pow(x, i);  // O(n^2)
    }
    /*calculo da derivada*/
    if (dpx != NULL) {
        for (int i = p.grau; i > 0; i--) {
            *dpx += p.p[i] * i * pow(x, i - 1);  // O(n^2)
        }
    }
}
