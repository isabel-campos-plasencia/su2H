
/* SU2 - Higgs   Header file */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>
#define DEBUG    
/*#define HISTER */  
#define twopi  6.28318
#define pi 3.14159
#define PUROGAUGE 

#define L    6

#define V (L*L*L*L)

#define maxit 5000    /* maximo numero de iteraciones por bin */
#define n_obs 3       /* numero de observables */

typedef struct
{
    float a0,a1,a2,a3;
} o4v;

/* random number generator */

#define NormRAN (1.0F/( (float) RAND_MAX+300.0F))
#define  RAN() ( (float) rand() * NormRAN )


#define FNORM   (2.3283063671E-10F)
#define RANDOM  ( (ira[ip++]=ira[ip1++]+ira[ip2++]) ^ira[ip3++] )
#define FRANDOM (FNORM*RANDOM)


#ifdef MAIN
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
#else
extern unsigned char ip,ip1,ip2,ip3;
extern unsigned ira[];
#endif

#define randmax NormRAN

#define randmax NormRAN
struct s_datos
{
    int itmax,        /*   Numero de medidas por bloque               */
    mfresh,           /* frecuencia de las medidas                    */
    nbin,            /* numero de bloque                              */
    itcut,           /* proximo bloque a calcular                     */
    flag,            /* conf de partida: 0(random), 1(fria),2(backup) */
    seed;            /* semilla random                                */
    float kappa1,         /* coupling a primeros vecinos               */
    kappa2,         /* coupling a segundos vecinos                */
    beta            /* coupling puro gauge   */     
#ifdef HISTER
        ,dkappa2,        /* variacion de kappa2 para la histeresis   */
         dbeta          /* variacion de beta para la histeresis      */
#endif
             ;   
};

#define NDATINT   6     /* numero de campos int   en s_datos */
#ifdef HISTER
#define NDATFLOAT 5     /* numero de campos float en s_datos */
#else
#define NDATFLOAT 3     /* numero de campos float en s_datos */
#endif

#define Normaenerp  ( (float) (1.0) / (6.0F * V) )    /* normalizaciones */
#define Normaener1  ( (float) (1.0) / (4.0F * V) )
#define Normaener2  ( (float) (1.0) / (6.0F * V) )

#define NOBS_HISTER 3
#define LPATH 100

FILE *Finput,*Foutput,*Fconfig;
o4v tirar;


