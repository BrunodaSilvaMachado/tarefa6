#include <stdio.h>
#include <linalg.h>

#define DEBUG 1
#define N 2

/**
*metodo de QuasiNewton para sistemas não lineares
*solução das questões 1a e 1d da seção 10.3
*/
//funções

//quetao 1a
double f1(double x[N])
{
	return 4*x[0]*x[0]  -20*x[0]+0.25*x[1]*x[1]+8;
}


double f2(double x[N])
{
	return 0.5*x[0]*x[1]*x[1] + 2*x[0]-5*x[1]+8;
}

//questao 1d
double f3(double x[N])
{
	return log((x[0]*x[0] + x[1]*x[1])/(2.0 * M_PI)) -sin(x[0]*x[1]);
}

double f4(double x[N])
{
	return exp(x[0] - x[1]) + cos(x[0]*x[1]);
}


/**
*@title: diff
*@descript: calculates the derivative numerical central 4 order
*@param f(): função vetorial
*@param x: vetor do numerico
*@param i: coordenada a ser derivada
*/

double diff(double f(),double x[N],int i)
{
    double h,a,b,c,d,dif;
    double tmp;
    h = 1e-6;
    tmp = x[i];
    x[i] = tmp + 2 * h;
    a = -f(x);
    x[i] = tmp + h;
    b = 8 * f(x);
    x[i] = tmp - h;
    c = -8 * f(x);
    x[i] = tmp - 2 * h;
    d = f(x);
    dif = (a + b + c + d) / (12 * h);
    x[i] = tmp;
    return dif;
}

//funçao a parte

void imprimeRaiz(double *r, int dim)
{
    int i;

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

double **multplicaMatrizEscalar(double num, double **matriz, int lin,int col)
{
    double **mat = (double**)malloc(lin*sizeof(double*));
	int i,j;
    for( i = 0; i < lin; i++)
        mat[i] = (double*)malloc(col*sizeof(double));

    for( i = 0; i < lin; i++)
    {
        for( j = 0; j < col; j++)
        {
            mat[i][j] = num * matriz[i][j];
        }
    }

    return mat;
}

double *metodoQuasiNewtow(double (*equacao[N])(), double x[N])
{
	double eps = 1e-6; //precição
	double **jf; // matriz jacobina aumentada
	double fx[N];
	double *y;
	int i,j,c = 0;

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
				jf[i][j] = diff(equacao[i],x,j);
			}

			jf[i][N] = fx[i];
		}

		#if DEBUG
		
		printf("X^(%d) =  ",c);
		for(i = 0; i < N;i++)
		{
			printf("%lf\t",x[i]);
		}
		puts("\n");
		c++;
		
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
	double (*equacao2a[N])() = {f1,f2};
	double x1[] = {0,0};

	double (*equacao2d[N])() = {f3,f4};
	double x2[] = {2,2};

	puts("set 10.3 questao 1-a");
	puts("Resultado final");
	imprimeRaiz(metodoQuasiNewtow(equacao2a,x1),N);

	puts("set 10.3 questao 1-d");
	puts("Resultado final");
	imprimeRaiz(metodoQuasiNewtow(equacao2d,x2),N);

	return 0;
}
