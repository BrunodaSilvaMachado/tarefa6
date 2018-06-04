#include <stdio.h>
#include <linalg.h>

#define DEBUG 1
#define N 3

/**
*metodo de newton para sistemas não lineares
*solução das questões 2a e 2d da seção 10.2
*/
//funções

//quetao 2a
double f1(double x[N])
{
	return 3*x[0] -cos(x[1]*x[2]) -0.5;
}


double f2(double x[N])
{
	return 4*x[0]*x[0] -625*x[1]*x[1] +2*x[1]-1;
}

double f3(double x[N])
{
	return exp(-x[0]*x[1]) + 20*x[2] + (10*M_PI -3)/3.0;
}

//questao 2d
double f4(double x[N])
{
	return 10*x[0] -2*x[1]*x[1] + x[1]-2*x[2]-5;
}


double f5(double x[N])
{
	return 8*x[1]*x[1]+4*x[2]*x[2]-9;
}

double f6(double x[N])
{
	return 8*x[1]*x[2]+4;
}

//derivadas 2a

double f1x1(double *x)
{
    return 3;
}

double f1x2(double *x)
{
    return -x[2]*sin(x[1]*x[2]);
}

double f1x3(double *x)
{
    return -x[1]*sin(x[1]*x[2]);
}

double f2x1(double *x)
{
    return 8*x[0];
}

double f2x2(double *x)
{
    return -1250*x[1] + 2;
}

double f2x3(double *x)
{
    return 0;
}

double f3x1(double *x)
{
    return -x[1]*exp(-x[0]*x[1]);
}

double f3x2(double *x)
{
    return -x[0]*exp(-x[0]*x[1]);
}

double f3x3(double *x)
{
    return 20;
}

//derivadas 2d

double f4x1(double *x)
{
    return 10;
}

double f4x2(double *x)
{
    return -4*x[1]+1;
}

double f4x3(double *x)
{
    return -2;
}

double f5x1(double *x)
{
    return 0;
}

double f5x2(double *x)
{
    return 16*x[1];
}

double f5x3(double *x)
{
    return 8*x[2];
}

double f6x1(double *x)
{
    return 0;
}

double f6x2(double *x)
{
    return 8*x[2];
}

double f6x3(double *x)
{
    return 8*x[1];
}


//funçao a parte

void imprimeRaiz(double *r, int dim)
{
    int i;

    puts("________");

    for(i = 0; i < dim; i++)
    {
        printf("x[%d] = %6.5lf\n",i + 1,r[i]);
    }

    puts("________");
}

double max(double * v, int n)
{
    int i;
    double m = v[0];

    for(i = 0; i < n;i++)
    {
        if(m < v[i])
            m = v[i];
    }

    return m;
}

double *metodoNewtow(double (*equacao[N])(), double (*jacobiano[N][N])(),double x0[N])
{
    double *x; //vetor tentativa
	double eps = 1e-6; //precição
	double **jf; // matriz jacobina aumentada
	double fx[N];
	double *y;
	int i,j;

	x = malloc(N * sizeof(double));

	memcpy(x,x0,N * sizeof(double));

    jf = (double**)malloc(N * sizeof(double*));

	for(i = 0; i < N; i++)
	{
		jf[i] = (double*)malloc((N + 1) * sizeof(double));
	}

	do
	{
		for(i = 0; i < N; i++)
		{
			fx[i] = -equacao[i](x);
			for(j = 0; j < N; j++)
			{
				jf[i][j] = (jacobiano[i][j])(x);
			}

			jf[i][N] = fx[i];
		}

		#if DEBUG
		puts("F(x)\n");
		imprimeRaiz(fx,N);

		puts("J(x)\n");
		imprimeMatriz(jf,N,N);
		#endif // DEBUG
		y = metodoGauss(jf,N,N + 1);

		for(i = 0; i < N; i++)
		{
			x[i] += y[i];
		}


	}while(max(y,N) > eps);

	return x;
}

int main(int argc, char **argv)
{

    double x[] = {0.1,0.1,-0.1};//vetor tentativa
	double (*equacao2a[N])() = {f1,f2,f3};
	double (*jacobiano2a[N][N])() = {{f1x1,f1x2,f1x3},{f2x1,f2x2,f2x3},{f3x1,f3x2,f3x3} };

	double (*equacao2d[N])() = {f4,f5,f6};
	double (*jacobiano2d[N][N])() = {{f4x1,f4x2,f4x3},{f5x1,f5x2,f5x3},{f6x1,f6x2,f6x3} };

	puts("set 10.2 questao 2-a");
	puts("Resultado final");
	imprimeRaiz(metodoNewtow(equacao2a,jacobiano2a,x),N);

	puts("set 10.2 questao 2-d");
	puts("Resultado final");
	imprimeRaiz(metodoNewtow(equacao2d,jacobiano2d,x),N);

	return 0;
}


