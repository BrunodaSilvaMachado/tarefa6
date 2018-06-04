/*
 * linalg.h
 *
 * Copyright 2018 Bruno da Silva Machado <brunosilvamachado@id.uff.br>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 */


#ifndef __LINALG_H__
#define __LINALG_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern double **lerMatriz(const char *, int *, int *);

extern void imprimeMatriz(double **, int, int );

extern double determinante(double **, int );

extern int triangularSuperior(double **, int, int );

extern int compararMatriz(double **, double **,int, int );

extern double *metodoJacobi(double **, int, double );

extern double *metodoGauss(double **, int, int  );

extern double **multiplicaMatriz(double **,int, int, double **, int, int );

extern double **somaMatriz(double **, double **,int, int );

extern double **diferencaMatriz(double **, double **,int, int );

extern double **matrizNula(int, int );

extern double **matrizInversa(double **, int, int );

extern double **matrizIdentidade(int, int );

extern double **matrizTransposta(double **, int, int );

extern double normalize(double *,int );

extern double *substituicaoRegressiva(double **, int );

#endif //LINALG_H__
