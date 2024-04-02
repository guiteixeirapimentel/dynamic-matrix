#pragma once
#include <math.h>
#include <string>
#include "Matriz.hpp"

Matriz CriaMatrizAleatoria(size_t numlinhas, size_t numColunas);

Matriz SubstSucessivas(const Matriz& L, const Matriz& b);

Matriz SubstRetroativas(const Matriz& U, const Matriz& b);

void DecompPALU(const Matriz& A, Matriz& POut, Matriz& LOut, Matriz& UOut);

void MostrarMatriz(const Matriz& M);

void LogMatrixToFile(const std::string& filename, const Matriz& M);

Matriz ResSistLinearPALU(const Matriz& A, const Matriz& b);

Matriz CalcInvMatriz(const Matriz& A);

void TrocaLinha(Matriz& m, size_t l1, size_t l2);

Matriz MatrizI(size_t n);
Matriz MatrizZeros(size_t n, size_t m);

size_t AchaIndicePivo(const Matriz& m, size_t nColuna, size_t nLinInicial);