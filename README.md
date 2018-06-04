# Sistemas não lineares

## RESOLUÇÃO DE SISTEMAS NÃO LINEARES
Uma equação que contenha uma expressão do tipo x^2, y^-2, x.y, √y/z, sen(x), exp(x+z), etc, é chamada não-linear em x, y, z, ..., porque ela não pode ser escrita como:
> 						ax + by + cz + ... = cte
que é uma equação linear em x, y, z, ...
Um sistema de *n* equações e *n* incógnitas x_1, x_2, ..., x_n é chamado de não-linear se uma ou mais equações é não-linear. Trazendo todos os termos diferentes de zero à esquerda de todas as equações, tem-se uma forma geral que pode ser usada para qualquer sistema não-linear.

![eq:1](imagens/01.png)

Em notação vetorial, o sistema linear acima pode ser escrito como: **F(x) = 0**, em que

![eq:2](imagens/02.png)

Um vetor X* = (x_1,x_2,...,x_n) que satisfaz **F(X\*) = 0** é denominado *raiz do sistema não linear*

# Metodos Utilizados

## Ponto fixo no R^n

No caso de um sistema não linear com N variaveis o procedimento se dá na seguinte forma, dado um conjunto de aproximações inicial x_0 para o vetor **x**, o método consiste em iterar sucessivamente a função dada sobre x_0, ou seja, constrói-se a sequência recursiva x_{n + 1} = f^{n + 1}(x_0) = f^{n}((f(x_0)) sendo cada x_{n} uma nova aproximação do ponto fixo **x** quando **F(x) = 0**, encerramos o loop e temos que **x** é a solução do sistema

## Método de Newton 

O método mais amplamente estudado e conhecido para resolver sistemas de equação não lineares é o Método de Newton.

O método se baseia em fazer uma escolha de uma aproximação inicial e atribuir em um vetor **x**, apos isso calcula-se **F(x)** e a matriz jacobiana **J(x)** associada a **F**. Depois resolvese o seguinte sistema linear

>	**J(x)y = -F(x)**

onde o vetor **y** é a icognita a ser determinada após isso atualizamos o valor de **x** definindo **x = x + y**, e por fim verificamos se a norma de **y** é menor que a tolerancia permitida, caso seja o programa se encerra e o valor **x** calculado é a aproximação da solução do sistema. 

A vantagem é que, sob certas condições sobre a aproximação inicial **x^(0)**, a função **F** e a matriz Jacobiana **J**, a sequência {x^(k)} produzida pelo método de Newton converge para a solução de **F(x) = 0** com taxa quadrática.

## Steepest Descents

O método Steepest Descents converge apenas linearmente para a solução, mas geralmente convergem mesmo para aproximações iniciais pobre. Como consequência, este método é usado para encontrar uma aproximação inicial com grande precisas para o metodo de Newton. O método de Steepest Descents determina um local mínimo para uma função de vvarias variável na forma g: R^n -> R. O método é valioso por ser independentemente da aplicação como do ponto de partida para resolver sistemas não-linear.

O metdo consiste nos seguintes passos

* Escolher uma aproximação inicial para o vetor **x**
* Calcular g onde g: R^n -> R
* Calcular o gradiente de g ,(z = grad g(**x^(k)**))
* obter a norma de g ,(z_0 = ||z||)
* calcular o vetor unitario ,(z = z/z_0)
* Determine a direção de **x^(0)** que diminui g, sugestão use o polinomio interpolador de Newton e apartir deste calcule o minimo
* Andar nesta direção e obtém x^(1), x = x -*a*z (onde *a* é o ponto de minimo do polinomio interpolado de Newton)
* Repetir até alcançar a precisão desejada


# Exercicios

## 1 metodo do Ponto fixo
![ex1](imagens/ex10.1.8.a.png)

``` 
x^(0) =  1.000000	1.000000	1.000000	

Resultado final
x[1] =  0.000000
x[2] =  2.000000
x[3] =  -1.000000
```

## 2 metodo de Newton
![ex2](imagens/ex10.2.2.a.png)
``` 
set 10.2 questao 2-a

X^(0) =  0.100000	0.100000	-0.100000	

F(x)

________
x[1] = 1.19995
x[2] = 7.01000
x[3] = -8.46203
________
J(x)

 3.000	-0.001	 0.001	
 0.800	-123.000	 0.000	
-0.099	-0.099	20.000	
________
X^(1) =  0.500106	0.045610	-0.521390	

F(x)

________
x[1] = -0.00060
x[2] = 1.20855
x[3] = -0.02163
________
J(x)

 3.000	-0.012	 0.001	
 4.001	-55.013	 0.000	
-0.045	-0.489	20.000	
________
Resultado final
________
x[1] = 0.49982
x[2] = 0.02362
x[3] = -0.52301
________
```
![ex3](imagens/ex10.2.2.d.png)
``` 
set 10.2 questao 2-d

X^(0) =  0.100000	0.100000	-0.100000	

F(x)

________
x[1] = 3.72000
x[2] = 8.88000
x[3] = -3.92000
________
J(x)

10.000	 0.600	-2.000	
 0.000	 1.600	-0.800	
 0.000	-0.800	 0.800	
________
X^(1) =  0.360000	6.300000	1.200000	

F(x)

________
x[1] = 76.88000
x[2] = -314.28000
x[3] = -64.48000
________
J(x)

10.000	-24.200	-2.000	
 0.000	100.800	 9.600	
 0.000	 9.600	50.400	
________
X^(2) =  0.524063	3.248634	0.501848	

F(x)

________
x[1] = 18.62167
x[2] = -76.43636
x[3] = -17.04255
________
J(x)

10.000	-11.995	-2.000	
 0.000	51.978	 4.015	
 0.000	 4.015	25.989	
________
X^(3) =  0.575809	1.811589	0.068083	

F(x)

________
x[1] = 4.13019
x[2] = -17.27338
x[3] = -4.98671
________
J(x)

10.000	-6.246	-2.000	
 0.000	28.985	 0.545	
 0.000	 0.545	14.493	
________
Resultado final
________
x[1] = 0.55598
x[2] = 1.22170
x[3] = -0.25383
________

```

## 3 metodo de Quasi-Newton
![ex4](imagens/ex10.3.1.a.png)
``` 
set 10.3 questao 1-a
X^(0) =  0.000000	0.000000	

F(x)

x[1] = -8.00000
x[2] = -8.00000
________
J(x)

-20.000	 0.000	
 2.000	-5.000	
________
X^(1) =  0.400000	1.760000	

F(x)

x[1] = -1.41440
x[2] = -0.61952
________
J(x)

-16.800	 0.880	
 3.549	-4.296	
________
X^(2) =  0.495894	1.983423	

F(x)

x[1] = -0.04926
x[2] = -0.05008
________
J(x)

-16.033	 0.992	
 3.967	-4.016	
________
X^(3) =  0.499988	1.999937	

F(x)

x[1] = -0.00014
x[2] = -0.00020
________
J(x)

-16.000	 1.000	
 4.000	-4.000	
________
X^(4) =  0.500000	2.000000	

F(x)

x[1] = -0.00000
x[2] = -0.00000
________
J(x)

-16.000	 1.000	
 4.000	-4.000	
________

Resultado final
x[1] = 0.50000
x[2] = 2.00000
________
```
![ex5](imagens/ex10.3.1.d.png)
``` 
set 10.3 questao 1-d

X^(0) =  2.000000	2.000000	

F(x)

x[1] = -0.99837
x[2] = -0.34636
________
J(x)

 1.807	 1.807	
 2.514	 0.514	
________

Resultado final
x[1] = 1.96868
x[2] = 1.47891
________

```

## 4 metodo Steepest Descents
![ex6](imagens/ex10.4.2.a.png)
``` 
set 10.4 questao 2-a

X^(0) =  0.242246	0.184829	0.371773	

Z(x)

-0.284493, -0.169658, -0.943546, 

X^(1) =  0.496634	0.333224	0.775835	

Z(x)

-0.508775, -0.296790, -0.808123, 

X^(2) =  0.913934	0.598953	0.848282	

Z(x)

-0.834600, -0.531459, -0.144894, 

X^(3) =  1.155587	0.883783	0.907899	

Z(x)

-0.638859, -0.753006, -0.157611, 

X^(4) =  1.029861	0.985397	0.935601	

Z(x)

0.766566, -0.619556, -0.168899, 

X^(5) =  1.043470	1.009093	0.909880	

Z(x)

-0.362668, -0.631427, 0.685399, 

X^(6) =  1.047267	1.031917	0.932791	

Z(x)

-0.116599, -0.700945, -0.703620, 

X^(7) =  1.048170	1.043737	0.920801	

Z(x)

-0.053567, -0.701033, 0.711114, 

X^(8) =  1.043938	1.055768	0.932300	

Z(x)

0.246445, -0.700605, -0.669640, 

Resultado final
1.043790, 1.062234, 0.925463, 
```
![ex6](imagens/ex10.4.2.d.png)
``` 
set 10.4 questao 2-d

X^(0) =  -0.027538	0.056913	0.381537	

Z(x)

0.255076, 0.086175, -0.963073, 

X^(1) =  0.024541	0.045191	0.878679	

Z(x)

-0.104158, 0.023444, -0.994284, 

Resultado final
-0.018456, 0.097891, 0.993464, 


```

para compilar os exercicios a seguir use:

>  gcc -L.\linalg\lib -Wl,-rpath=.\linalg\lib -llinalg -o main .\set10.1\main.c -I.\linalg\include

# menções

@thadeupenna

# Referencias

wwwp.fc.unesp.br\Numerico\SNLinear
