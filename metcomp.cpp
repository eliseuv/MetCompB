#include "MyLib.h"

//===========================================================================================================
// MetCompA
//===========================================================================================================

/*
Derivadas da funcao f(x)

Argumentos:
	*f = ponteiro da funcao
	x = ponto
	dx = elemento infinitesimal
*/

long double deriv_dir(long double (*f)(long double), long double x, long double dx){

return (f(x+dx)-f(x))/dx;

}

long double deriv_cent(long double (*f)(long double), long double x, long double dx){

return (f(x+dx)-f(x-dx))/(2.0*dx);

}


long double deriv_seg(long double (*f)(long double), long double x, long double dx){

return (f(x+dx)+f(x-dx)-2.0*f(x))/(dx*dx);

}

/*
Integrais da funcao f(x)

Argumentos:
	*f = ponteiro da funcao
	a = inicio do intervalo
	b = final do intervalo
	n = numero elementos infinitesimais
*/

long double int_esq(long double (*f)(long double), long double a, long double b, int n){

long double dx, ie = 0;
int i;

dx=(b-a)/n;

for (i = 0; i < n; i++){
	ie += f(a+i*dx);
	}
ie *= dx;

return ie;

}

long double int_cent(long double (*f)(long double), long double a, long double b, int n){

long double dx, ic = 0;
int i;

dx=(b-a)/n;

for (i = 0; i < n; i++){
	ic += f(a+(i+0.5)*dx);
	}
ic *= dx;

return ic;

}

long double int_trap(long double (*f)(long double), long double a, long double b, int n){

long double dx, it = 0;
int i;

dx=(b-a)/n;

for (i = 0; i < n; i++){
	it += f(a+(i+1)*dx)+f(a+i*dx);
	}
it *= (dx)/2;;

return it;

}

long double int_simp(long double (*f)(long double), long double a, long double b, int n){

long double dx, is = 0;
int i;

dx=(b-a)/n;

for (i = 0; i < n; i++){
	is += f(a+i*dx)+4*f(a+(i+0.5)*dx)+f(a+(i+1)*dx);
	}
is *= dx/6;

return is;

}

/*
Zeros da funcao f(x)

Argumentos:
	*f = ponteiro da funcao
	*fp = ponteiro da funcao derivada
	x, x1, x2 = chutes iniciais
	dx = elemento infinitesimal
	d = erro para parada
	max = maximo do contador
*/

long double zero_iteracao(long double (*f)(long double), long double x, long double d, int max){

int cont = 0;

while ((cont<max)&&(abs(f(x))>d)){
	cont++;
	x=x-f(x); 
	}

return x;

}

long double newton_raphson(long double (*f)(long double), long double (*fp)(long double), long double x, long double d, int max){

int cont = 0;

while ((cont<max)&&(abs(f(x))>d)){
	cont++;
	x=x-f(x)/fp(x); 
	}

return x;

}

long double newton_raphson(long double (*f)(long double), long double x, long double d, long double dx, int max){

int cont = 0;

while ((cont<max)&&(abs(f(x))>d)){
	cont++;
	x=x-f(x)/deriv_cent(f,x,dx); 
	}

return x;

}

long double bisseccao (long double (*f)(long double), long double x1, long double x2, long double d, int max){

int cont = 0;
long double xm=(x1+x2)/2;

while ((cont<max)&&(abs(f(xm))>d)){
	cont++;
	xm=(x1+x2)/2;
	if (f(x1)*f(xm)>0){
		x1=xm;
	}
	else{
		x2=xm;
		}
	}
return xm;
}

// Interpolacao

/*
Metodo de Neville

Argumentos:
	n = numero de pontos
	X[] = array de x
	Y[] = array de y
	x = ponto

*/

long double neville(int n, long double X[], long double Y[], long double x){

long double P[n];
int i, j;

for (i=0;i<n;i++){
	P[i]=Y[i];
	}
for (j=n-1;j>0;j--){
	for (i=0;i<j;i++){
		P[i]=((x-X[i+n-j])*P[i]+(X[i]-x)*P[i+1])/(X[i]-X[i+n-j]);
		}
	}

return P[0];

}

/*
Spline cubico

Argumentos:
	N = numero de pontos
	x[] = array de x
	y[] = array de y
	X = ponto

*/

long double spline(int N, long double x[], long double y[], long double X){

const int n = N;
int i, l;
long double a[n-2], b[n-2], c[n-2], f[n-2], alpha[n-2], beta[n-2], gamma[n-1], y2[n];

alpha[n-2] = 0;
beta[n-2] = 0;
y2[0] = 0;
y2[n-1] = 0;

for (i = 1; i< n-1; i++){
	a[i]=(x[i]-x[i-1])/6;
	b[i]=(x[i+1]-x[i-1])/3;
	c[i]=(x[i+1]-x[i])/6;
	f[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
	}

for (i = n-2; i > 0; i--){
	gamma[i]=-1/(c[i]*alpha[i]+b[i]);
	alpha[i-1]=gamma[i]*a[i];
	beta[i-1]=gamma[i]*(c[i]*beta[i]-f[i]);
	}

for (i=0; i<n-2; i++){
	y2[i+1]=y2[i]*alpha[i]+beta[i];
}

for (i=0;i<n-1;i++){
if ((X>x[i])&&(X<x[i+1])){
	l=i;
	break;
	}
}

return (y[l]*((x[l+1]-X)/(x[l+1]-x[l]))+y[l+1]*((X-x[l])/(x[l+1]-x[l]))+y2[l]*(((x[l+1]-X)/(x[l+1]-x[l]))*(pow(((x[l+1]-X)/(x[l+1]-x[l])),2)-1)*pow(x[l+1]-x[l],2))/6+y2[l+1]*(((X-x[l])/(x[l+1]-x[l]))*(pow(((X-x[l])/(x[l+1]-x[l])),2)-1)*pow(x[i+1]-x[i],2))/6);

}

// Fit

/*
Minimos quadrados

y(x) = ax + b

Argumentos:
	n = numero de pontos
	x[] = array de x
	y[] = array de y
	*a = ponteiro de a
	*b = ponteiro de b
	k[] = array com coeficientes
	X = ponto

*/

void least_squares(int n, long double x[], long double y[], long double *a, long double *b){

	long double Sx=0, Sy=0, Sxx=0, Sxy=0;
	int i;

	for (i=0; i<n; i++){
		Sx += x[i];
		Sy += y[i];
		Sxx += x[i]*x[i];
		Sxy += x[i]*y[i];
		}

	if (n*Sxx-Sx*Sx != 0){
		*a = (n*Sxy-Sx*Sy)/(n*Sxx-Sx*Sx);
		*b = (Sy*Sxx-Sx*Sxy)/(n*Sxx-Sx*Sx);
		}
	else{
		cout << "Sistema impossivel!";
		}
	
}

void least_squares(int n, long double x[], long double y[], long double k[]){

	long double Sx=0, Sy=0, Sxx=0, Sxy=0, v[2];
	int i;

	for (i=0; i<n; i++){
		Sx += x[i];
		Sy += y[i];
		Sxx += x[i]*x[i];
		Sxy += x[i]*y[i];
		}

	if (n*Sxx-Sx*Sx != 0){
		k[0] = (n*Sxy-Sx*Sy)/(n*Sxx-Sx*Sx);
		k[1] = (Sy*Sxx-Sx*Sxy)/(n*Sxx-Sx*Sx);
		}
	else{
		cout << "Sistema impossivel!";
		}
	
}

long double least_squares(int n, long double x[], long double y[], long double X){

	long double Sx=0, Sy=0, Sxx=0, Sxy=0;
	int i;

	for (i=0; i<n; i++) {
		Sx += x[i];
		Sy += y[i];
		Sxx += x[i]*x[i];
		Sxy += x[i]*y[i];
		}
	
	if (n*Sxx-Sx*Sx != 0){
		return ((n*Sxy-Sx*Sy)*X+(Sy*Sxx-Sx*Sxy))/(n*Sxx-Sx*Sx);
		}
	else{
		cout << "Sistema impossivel!";
		return 0;
		}
}

//===========================================================================================================
// MetCompB
//===========================================================================================================

// Operadores extrapolados

/*
Extrapola y'(x) para dx=0

Argumentos:
	f = f(x)
	x = ponto para f'(x)
	dx = menor dx utilizado
	n = numero de pontos para extrapolacao
	
*/

long double derivex (long double (*f)(long double), long double x, long double dx, int n){
	
	long double X[n], Y[n], P[n], d;
	int i, j;
	
	for (i=0; i<n; i++){
		d = (i+1)*dx;
		X[i] = d;
		Y[i] = (f(x+d)-f(x-d))/(2*d);
	}

	for (i=0;i<n;i++){
		P[i]=Y[i];
		}
	for (j=n-1;j>0;j--){
		for (i=0;i<j;i++){
			P[i]=((-X[i+n-j])*P[i]+(X[i])*P[i+1])/(X[i]-X[i+n-j]);
			}
		}

	return P[0];
	
}

/*
Extrapola int_a^b f(x) dx para dx=0

Argumentos:
	f = f(x)
	a = ponto inicial
	b = ponto final
	n = maior numero de intervalos utilizados
	m = numero de pontos para extrapolacao
	
*/

long double intsimpex (long double (*f)(long double), long double a, long double b, int n, int m){

long double dx, is = 0;
long double X[m], Y[m], P[m];
int i, j, k = n;

for (j=1; j<=m; j++){
	dx=(b-a)/k;
	for (i = 0; i < k; i++){
		is += f(a+i*dx)+4*f(a+(i+0.5)*dx)+f(a+(i+1)*dx);
		}
	is *= dx/6;
	X[j] = dx;
	Y[j] = is;
	is = 0;
	if (k%2 == 0){
		k = k/2;
	}
	else{
		k = (k+1)/2;
	}
}

return neville(m,X,Y,0);

}

// EDO

/*
Método de Euler para EDO homogênea de primeira ordem

Resolve para f(p) em 
	
	f'(x) = g(x)

Argumentos:
	g = g(x)
	p = ponto para f(p)
	n = numero de passos
	a = ponto conhecido
	fa = f(a)
	
*/

long double eulerx (long double (*g)(long double), long double p, int n, long double a, long double fa){
	
	long double dx = (p-a)/n;
	int i;
	
	for (i=0; i<=n; i++){
		fa = fa + g(a)*dx;
		a = a + dx;
	}
	
	return fa;
	
}

/*
Método de Euler para EDO autônoma de primeira ordem

Resolve para f(p) em 
	
	f'(x) = g(f(x))

Argumentos:
	g = g(f(x))
	p = ponto para f(p)
	n = numero de passos
	a = ponto conhecido
	fa = f(a)
	
*/

long double eulery (long double (*g)(long double), long double p, int n, long double a, long double fa){
	
	long double dx = (p-a)/n;
	int i;
	
	for (i=0; i<=n; i++){
		fa = fa + g(fa)*dx;
	}
	
	return fa;
	
}

/*
Método de Euler para EDO de primeira ordem

Resolve para y(p) em 
	
	y'(x) = g(x,y(x))

Argumentos:
	g = g(x,y(x))
	p = ponto para y(p)
	n = numero de passos
	a = ponto conhecido
	ya = y(a)
	
*/

long double euler (long double (*g)(long double, long double), long double p, int n, long double a, long double ya){
	
	long double dx = (p-a)/n;
	int i;
	
	for (i=0; i<=n; i++){
		ya += g(a,ya)*dx;
		a += dx;
	}
	
	return ya;

}

/*
Método Runge-Kutta de segunda ordem para EDO de primeira ordem

Resolve para y(p) em 
	
	y'(x) = g(x,y(x))

Argumentos:
	g = g(x,y(x))
	p = ponto para y(p)
	n = numero de passos
	a = ponto conhecido
	ya = y(a)
	
*/

long double rungekutta2 (long double (*g)(long double, long double), long double p, int n, long double a, long double ya){
	
	long double dx = (p-a)/n;
	int i;
	
	for (i=0; i<=n; i++){
		ya += dx*g(a+(dx/2),ya+(dx/2)*g(a,ya));
		a += dx;
	}
	
	return ya;
}

/*
Método Runge-Kutta de quarta ordem para EDO de primeira ordem

Resolve para y(p) em 
	
	y'(x) = g(x,y(x))

Argumentos:
	g = g(x,y(x))
	p = ponto para y(p)
	n = numero de passos
	a = ponto conhecido
	ya = y(a)
	
*/

long double rungekutta4 (long double (*g)(long double, long double), long double p, int n, long double a, long double ya){
	
	long double k1, k2, k3, k4, dx = (p-a)/n;
	int i;
	
	for (i=0; i<=n; i++){
		
		k1 = dx*g(a,ya);
		k2 = dx*g(a+dx/2,ya+k1/2);
		k3 = dx*g(a+dx/2,ya+k2/2);
		k4 = dx*g(a+dx,ya+k3);
		
		ya += k1/6 + k2/3 + k3/3 + k4/6;
		
		a += dx;
	}
	
	return ya;

}

/*
Método Ponto médio modificado para integrar EDO de primeira ordem

Resolve para y(p) em 
	
	y'(x) = g(x,y(x))

Argumentos:
	g = g(x,y(x))
	p = ponto para y(p)
	n = numero de passos
	a = ponto conhecido
	ya = y(a)
	
*/

long double ptmedio (long double (*g)(long double, long double), long double p, int n, long double a, long double ya){
	
	long double dx=(p-a)/n, z[3], ym;
	int i;
	
	if(n%2==1){n++;}
	
	z[0] = ya;
	z[1] = z[0] + dx*g(a,z[0]);
	
	for (i=0; i<(n-2)/2; i++){
		z[(i+2)%3] = z[i%3] + 2*dx*g(a+(i+2)*dx,z[(i+1)%3]);
	}
	
	ym = 0.5*(z[((n-2)/2+2)%3] + z[((n-2)/2+1)%3] + dx*g(a+n*dx,z[((n-2)/2+2)%3]));
	
	cout << ym << endl;
	
	for (i=(n-2)/2+1; i<n; i++){
		z[(i+2)%3] = z[i%3] + 2*dx*g(a+(i+2)*dx,z[(i+1)%3]);
	}
	
	ya = 0.5*(z[(n+1)%3] + z[n%3] + dx*g(a+n*dx,z[(n+1)%3]));
	
	cout << ya << endl;
	
	return (4*ya-ym)/3;

}

//===========================================================================================================
// Outros
//===========================================================================================================

// File

/*
Escrever

Argumentos:
	f1 = funcao 1
	f2 = funcao 2
	a = inicio do intervalo
	b = final do intervalo
	delta = incremento

*/

void write (long double (*f1)(long double), long double (*f2)(long double) , int a, int b, float delta){

	int i;

	ofstream txt("data.txt", ios::out);

	for (i = a; i < b; i++){
		txt << f1(i*delta) << "\t" << f2(i*delta) << endl;
		cout << f1(i*delta) << "\t" << f2(i*delta) << endl;
		}
}

// Gauss-Jordan elimination

void matrix_build (int m, int n, long double **A){
	
	int i;

	A = new long double *[m];
		for(i = 0; i < m; i++)
    		A[i] = new long double[n];
	}

void gauss_jordan (int m, int n, int o, RealMatrix A, RealMatrix B){

	int i,j,k, rank = min(m,n);
	long double v;	// Aux
	
	// Triangularization
	for (i=0; i<n-1; i++){	// Cols loop 1
		// Partial pivoting
		for (j=i+1; j<n; j++){	// Lines loop
			if (abs(A[j][i])>abs(A[i][i])){	// Pivot test
				for (k=0; k<n; k++){	// Cols loop
					v = A[i][k];		//
					A[i][k] = A[j][k];	// Bubble sort
					A[j][k] = v;		//
					}
				for (k=0; k<m; k++){	// Cols loop
					v = B[i][k];		//
					B[i][k] = B[j][k];	// Bubble sort
					B[j][k] = v;		//
					}
				}
			}
		// 
		for (j=i+1; j<n; j++){	// Lines loop
			v = A[j][i]; 	// Store Aji element
			for (k=0; k<n; k++){	// Cols loop A
				A[j][k]-=(v/A[i][i])*A[i][k];
				}
			for (k=0; k<m; k++){	// Cols loop B
				B[j][k]-=(v/A[i][i])*B[i][k];
				}
			}
		}
	// Diagonal
	for (i=n-1; i>0; i--){	// Cols loop 1
		for (j=i-1; j>-1; j--){	// Lines loop
			v = A[j][i]; 	// Store Aji element
			for (k=0; k<n; k++){	// Cols loop A
				A[j][k]-=(v/A[i][i])*A[i][k];
				}
			for (k=0; k<m; k++){	// Cols loop B
				B[j][k]-=(v/A[i][i])*B[i][k];
				}
			}
		}
	// Normalization
	for (i=0; i<n; i++){	// Diagonal loop
		v = A[i][i];	// Store diagonal element
		for (j=0; j<n; j++){	// Cols loop
			A[i][j] /= v;	
			}
		for (j=0; j<m; j++){	// Cols loop
			B[i][j] /= v;	
			}
		}
}

// Quantum algorithms

/*
	Quantum Phase Estimation Routine

	Argumentos:
		p = phase
		t = number of qubits in the first register
		l = number of estimatives
		ev[] = array with estimatives values
		ep[] = array with estimatives probs

*/

void QPE (long double p, int t, int l, long double ev[], long double ep[]){
	
	const long double pi = 4.0*atan(1.0);
	int i,j;
	
	int n = static_cast<int> (pow(2,t)); // n = 2^t
	vector<complex<long double> > r(n,0);	// First register
	long double P = n*p, A = 2*pi/n;
	
	// Apply circuit
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			r[i] += (cos(A*j*(P-i)),sin(A*j*(P-i)));
		}
		r[i] /= n;
	}
	
	// Bubble-sort
	int index [l];
	complex<long double> k;
	for (i = 0; i < l; i++){
		for (j = i; j < n; j++){
			if (abs(r[j])>abs(r[i])){
				k = r[i];
				r[i] = r[j];
				r[j] = k;
				index[i] = j;
			}
		}
	}
	
	// Write results
	for (i = 0; i < l; i++){
		ev[i] = static_cast<long double>(index[i])/n;
		ep[i] = abs(r[i]);
	}
	
}

void QPE (long double p, int m, long double e, int l, long double ev[], long double ep[]){
	
	const long double pi = 4.0*atan(1.0);
	int i,j;
	
	int t = m + static_cast<int>(ceil(log2(2+1/(2*e))));
	
	int n = static_cast<int> (pow(2,t)); // n = 2^t
	vector<complex<long double> > r(n,0);	// First register
	long double P = n*p, A = 2*pi/n;
	
	// Apply circuit
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			r[i] += (cos(A*j*(P-i)),sin(A*j*(P-i)));
		}
		r[i] /= n;
	}
	
	// Bubble-sort
	int index [l];
	complex<long double> k;
	for (i = 0; i < l; i++){
		for (j = i; j < n; j++){
			if (abs(r[j])>abs(r[i])){
				k = r[i];
				r[i] = r[j];
				r[j] = k;
				index[i] = j;
			}
		}
	}
	
	// Write results
	for (i = 0; i < l; i++){
		ev[i] = static_cast<long double>(index[i])/n;
		ep[i] = abs(r[i]);
	}
	
}