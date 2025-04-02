#include "SistemaLinear.h"
#include <stdlib.h>
#include <stdio.h>
#include "utils.h"

// #include "likwid.h"

int main(){
    int tamanho = 0;
    scanf("%d", &tamanho);
    s_linear *sistema = cria_s_linear(tamanho);
    // LIKWID_MARKER_INIT;
    // LIKWID_MARKER_START ("Teste_LIKWID");
    // //calculo
    // LIKWID_MARKER_STOP ("Teste_LIKWID");
    // LIKWID_MARKER_CLOSE;
    return 0;
}