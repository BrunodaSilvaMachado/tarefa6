#include <stdio.h>
#include <linalg.h>

#define DEBUG 0

#define N 3
#define MAX 1000
#define TOL 0.05

/**
*metodo de Steepest Descent para sistemas não lineares
*solução das questões 2a e 2d da seção 10.4
*/

double f1(double x[N])
{
    return 15*x[0] + x[1]*x[1] -4*x[2] - 13;
}

double f2(double x[N])
{
    return x[0]*x[0]+10*x[1]-x[2] -11;
}

double f3(double x[N])
{
    return x[1]*x[1]*x[1] -25*x[2] + 22;
}

double f4(double x[N])
{
    return x[0] + cos(x[0]*x[1]*x[2]) -1;
}

double f5(double x[N])
{
    return pow(1-x[0],0.25) + x[1] +0.05*x[2]*x[2] -0.15*x[2] -1;
}

double f6(double x[N])
{
    return -x[0]*x[0] -0.1*x[1]*x[1] + 0.01*x[1] + x[2] -1;
}

/**
*@descript gn:Soma dos modulo ao quadrado das funções vetoriais
*
*/
double gn(double (*f[N])(), double x[N])
{
    return pow(f[0](x),2) + pow(f[1](x),2) + pow(f[2](x),2);
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

double **multplicaMatrizEscalar(double num, double **matriz, int lin,int col)
{
    double **mat = (double**)malloc(lin*sizeof(double*));
    for(int i = 0; i < lin; i++)
        mat[i] = (double*)malloc(col*sizeof(double));

    for(int i = 0; i < lin; i++)
    {
        for(int j = 0; j < col; j++)
        {
            mat[i][j] = num * matriz[i][j];
        }
    }

    return mat;
}

double *multplicaVetorEscalar(double num, double *vetor, int col)
{
    double *vet = (double*)malloc(col * sizeof(double));
    for(; col > 0; col--)
    {
        vet[col - 1] = num * vetor[col -1];
    }

    return vet;
}


double *grad(double (*equacao[N])(), double x[N])
{
    double **fx;
    double **result;
    double **jf;
    int i,j;

    fx = (double**)malloc(N*sizeof(double*));
    jf = (double**)malloc(N*sizeof(double*));

    for(i = 0; i < N; i++)
    {
        fx[i] = (double*)malloc(1*sizeof(double));
        jf[i] = (double*)malloc(N*sizeof(double));
    }

    for(i = 0; i < N; i++)
    {
        fx[i][0] = equacao[i](x);
        for(j = 0; j < N; j++)
        {
            jf[i][j] = diff(equacao[i],x,j);
        }
    }
    result = multplicaMatrizEscalar(2,multiplicaMatriz(matrizTransposta(jf,N,N),N,N,fx,N,1),N,1);
    result = matrizTransposta(result,N,1);

#if DEBUG
    imprimeMatriz(jf,N,N);
    imprimeMatriz(result,1,N);
#endif // DEBUG

    return *result;

}

void imprimeVetor(double *v, int dim)
{
    for(int i = 0; i < dim; i++)
    {
        printf("%lf, ",v[i]);
    }
    puts("\n");
}

double *subtraiVetor(double *a, double *b, int dim)
{
    double *sub = (double*)malloc(dim * sizeof(double));

    for(int i = 0; i < dim; i++)
    {
        sub[i] = a[i] - b[i];
    }

    return sub;
}

double *steepestDescent(double (*f[N])(), double x0[N])
{
    double *z; //vetor gradiente normatizado
    double z0; //norma do gradiente
    double *g; //vetor gradiente
    double *x; //vetor tentativa
    /**polinomio interpolador de newton*/
    double g1,g2,g3;
    double a,a0,a1,a2,a3;
    double h1,h2,h3;

    int k = 0;//contador
    x = (double*)malloc(N * sizeof(double));

    memcpy(x,x0,N * sizeof(double));

    while(k < MAX)
    {
        #if DEBUG
        puts("step 3");
        #endif // DEBUG
        g1 = gn(f,x);
        g = grad(f,x);
        z0 = normalize(g,N);

        #if DEBUG
        puts("step 4");
        #endif // DEBUG
        if(z == 0)
            break;

        #if DEBUG
        puts("step 5");
        #endif // DEBUG
        z = multplicaVetorEscalar(1.0/z0,g,N);
        #if DEBUG
        puts("vet z");

        imprimeVetor(z,N);
        #endif // DEBUG
        a1 = 0;
        a3 = 1;
        g3 = gn(f,subtraiVetor(x,multplicaVetorEscalar(a1,z,N),N));

        #if DEBUG
        puts("step 6");
        puts("vet6");
        imprimeVetor(x,N);
        #endif // DEBUG
        while(g3 >= g1)
        {
            #if DEBUG
            puts("step 7");
            #endif // DEBUG
            a3 = a3/2.0;
            g3 = gn(f,subtraiVetor(x,multplicaVetorEscalar(a3,z,N),N));

            #if DEBUG
            puts("step 8");
            #endif // DEBUG

            if(a3 < TOL/2)
            {
                break;
            }
        }

        #if DEBUG
        puts("step 9");
        #endif // DEBUG

        a2 = a3/2;
        g2 = gn(f,subtraiVetor(x,multplicaVetorEscalar(a2,z,N),N));
        #if DEBUG
        puts("step 10");
        #endif // DEBUG

        h1 = (g2 - g1)/(double)a2;
        h2 = (g3 - g2)/(double)(a3 - a2);
        h3 = (h2 - h1)/(double)a3;

        #if DEBUG
        puts("step 11");
        #endif // DEBUG
        a0 = 0.5*(a2 - h1/(double)h3);

        #if DEBUG
        puts("step 12");
        #endif // DEBUG
        a = (a0 < a3)?a0:a3;

        #if DEBUG
        printf("a = %lf\n",a);

        puts("step 13");

        #endif // DEBUG

        for(int i = 0; i < N; i++)
        {
            x[i] = x[i] - a*z[i];
        }

        #if DEBUG
        puts("vet13     x");
        imprimeVetor(x,N);
        puts("vet13     z");
        imprimeVetor(z,N);

        puts("step 14");

        #endif // DEBUG

        if(fabs(gn(f,x) - g1) < TOL)
            break;

        #if DEBUG
        puts("step 15");
        #endif // DEBUG
        k++;
    }

    if(k == MAX)
        puts("Maximo de interações exedido");

    return x;

}

int main()
{
    double x1[N] = {0.1,0.1,-0.1}; //vetor tentativa
    double (*equacao2a[N])(double *) = {f1,f2,f3}; //vetror de funções
    double (*equacao2d[N])(double *) = {f4,f5,f6}; //vetror de funções

    puts("set 10.4 questao 2-a");
	puts("Resultado final");
	imprimeVetor(steepestDescent(equacao2a,x1),N);

	puts("set 10.4 questao 2-d");
	puts("Resultado final");
	imprimeVetor(steepestDescent(equacao2d,x1),N);

    return 0;
}
