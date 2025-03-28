#include <float.h>
#include <math.h>
#include <stdio.h>

#include "ZeroFuncao.h"
#include "utils.h"

int main() {
    real_t a, b;
    Polinomio pol;
    printf("grau do polinomio:\n");
    scanf("%d", &pol.grau);
    pol.p = (real_t*)malloc((pol.grau + 1) * sizeof(int));
    printf("coeficientes do polinomio:\n");
    for (int i = pol.grau; i >= 0; --i) {
        scanf("%lf", &pol.p[i]);  // le os coeficientes do polinomio
    }
    printf("Intervalo das raizes:\n");
    scanf("%lf %lf", &a, &b);  // intervalo onde está uma das raizes.

    for (int i = 1; i < 4; i++) {
        int it = 0;
        real_t raiz = 0.0;
        real_t erro = bisseccao(pol, a, b, i, &it, &raiz);
        printf("bissec <raiz: %.15lf> <criterio(%d): %lf> <it: %d>\n", raiz, i, erro, it);
    }

    for (int i = 1; i < 4; i++) {
      int it = 0;
      real_t raiz = 0.0;
      real_t erro = newtonRaphson(pol, a, i, &it, &raiz);
      printf("newton <raiz: %.15lf> <criterio(%d): %lf> <it: %d>\n", raiz, i, erro, it);
  }
    return 0;
}
