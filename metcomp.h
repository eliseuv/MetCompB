#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <complex>
#include <vector>

using namespace std;

// Types

// Real matrix type
typedef vector<vector<long double> > RealMatrix;
// Complex matrix type
typedef vector<vector<complex<long double> > > ComplexMatrix;

// Functions

//===========================================================================================================
// MetCompA
//===========================================================================================================

// Derivadas
long double deriv_dir		(long double (*)(long double), long double, long double);
long double deriv_cent		(long double (*)(long double), long double, long double);
long double deriv_seg		(long double (*)(long double), long double, long double);

// Integrais
long double int_esq			(long double (*)(long double), long double, long double, int);
long double int_cent		(long double (*)(long double), long double, long double, int);
long double int_trap		(long double (*)(long double), long double, long double, int);
long double int_simp		(long double (*)(long double), long double, long double, int);

// Zeros
long double zero_iteracao	(long double (*)(long double), long double, long double, int);
long double newton_raphson	(long double (*)(long double), long double (*)(long double), long double, long double, int);
long double newton_raphson	(long double (*)(long double), long double, long double, long double, int);
long double bisseccao		(long double (*)(long double), long double, long double, long double, int);

// Interpolacao
long double neville			(int, long double [], long double [], long double);
long double spline			(int, long double [], long double [], long double);

// Fit
void least_squares			(int, long double [], long double [], long double * , long double *);
void least_squares			(int, long double [], long double [], long double []);
long double least_squares	(int, long double [], long double [], long double);

//===========================================================================================================
// MetCompB
//===========================================================================================================

// Operadores extrapolados
long double derivex			(long double (*)(long double), long double, long double, int);
long double intsimpex		(long double (*)(long double), long double, long double, int, int);

// EDO
// Euler
long double euler 			(long double (*)(long double, long double), long double, int, long double, long double);
long double eulerx 			(long double (*)(long double), long double, int, long double, long double);
long double eulery 			(long double (*)(long double), long double, int, long double, long double);

// Runge-Kutta
long double rungekutta2		(long double (*)(long double, long double), long double, int, long double, long double);
long double rungekutta4		(long double (*)(long double, long double), long double, int, long double, long double);

// Ponto medio modificado
long double ptmedio			(long double (*)(long double, long double), long double, int, long double, long double);

//===========================================================================================================
// Outros
//===========================================================================================================

// Escrever
void write					(long double (*)(long double), long double (*)(long double), int, int, float);

// Matrix
void matrix_build			(int, int, long double **);
void gauss_jordan			(int, int, int, RealMatrix, RealMatrix);

// Quantum
void QPE (long double, int, int, long double [], long double []);
void QPE (long double, int, long double, int, long double [], long double []);