#include <stdio.h>
#include <math.h>

/**
* Ponto fixo.
* solução da questão 8-d
*/

#define N 3

double f1(double x[N])
{
	return x[0] * x[0] + 2 * x[1] * x[1] - x[1] - 2 * x[2];
}


double f2(double x[N])
{
	return x[0] * x[0] - 8 * x[1] * x[1] + 10 * x[2];
}


double f3(double x[N])
{
	return x[0] * x[0]/(7.0 * x[1] * x[2]) -1;
}

double normalize(double *v,int dim)
{
	int i;
	double sum = 0;
	for(i = 0; i < dim; i++)
	{
		sum += v[i] * v[i];
	}

	sum = sqrt(sum);

	return sum;
}

int main(int argc, char **argv)
{
	double xa[N] = {1,1,1};
	double eps = 1e-5,tol,norm,norma;
	double (*equacao[N])() = {f1,f2,f3};

	int i,c = 0;

	do
	{
		printf("x^(%d) =  ",c);
		for(i = 0; i < N;i++)
		{
			printf("%lf\t",xa[i]);
		}
		puts("\n");
		c++;
		
		for(i = 0; i < N;i++)
		{
			norma = normalize(xa,N);

			xa[i] = equacao[i](xa);

			norm = normalize(xa,N);
		}
		
		tol = fabs(norm - norma)/norm;

	}while(tol > eps);

	puts("Resultado final");
	for(i = 0; i < N;i++)
	{
		printf("x[%d] =  %lf\n",i + 1,xa[i]);
	}

}
