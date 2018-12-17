 /* SU2H_INI.C */

#include "su2h.h"
#include "su2h_arit.h"
#include <math.h>

extern int x_p[],y_p[],z_p[],t_p[],
           x_m[],y_m[],z_m[],t_m[];

extern int jp[],jm[],neigh_p[],neigh_n[];


extern o4v f[],U[4][V];     
 
extern struct s_datos datos;

extern float obs_dat[n_obs][maxit],resultados[NOBS_HISTER];


extern float kappa1,kappa2,beta;

char dir[LPATH];


void Init_Rand(int semilla)     /* Se utiliza el generador estandar de C  */
{                               /* para inicializar otro generador random */
    int i;                      /* mas seguro. Para aumentar mas la segu- */
                                /* ridad se eliminan los primeros numeros */
    srand((unsigned)semilla);
    for (i=0;i<111;i++)         /* Probablemente no es necesario */
        rand();

    ip=128;    
    ip1=ip-24;    
    ip2=ip-55;    
    ip3=ip-61;
    
    for (i=0; i<256; i++)
        ira[i] = (unsigned) rand()+ (unsigned) rand(); 

    for (i=0;i<1111;i++)        /* Probablemente tampoco es necesario */
        RANDOM;
}


void Direccionamientos(void)
{
    int i;

    for (i=0;i<L;i++)
    {
        x_p[i]= 1;
        x_m[i]=-1;
        y_p[i]= L;
        y_m[i]=-L;
        z_p[i]= L*L;
        z_m[i]=-L*L;
        t_p[i]= L*L*L;
        t_m[i]=-L*L*L;
    }
    x_m[0]= L-1;
    y_m[0]=(L-1)*L;
    z_m[0]=(L-1)*L*L;
    t_m[0]=(L-1)*L*L*L;
    x_p[L-1]=-x_m[0];
    y_p[L-1]=-y_m[0];
    z_p[L-1]=-z_m[0];
    t_p[L-1]=-t_m[0];
} 

void Neigh(int j)
{
    

    jp[0] = j + neigh_p[0];        /* CALCULO DE LOS SITES VECINOS */
    jm[0] = j + neigh_n[0];
    jp[1] = j + neigh_p[1];
    jm[1] = j + neigh_n[1];        /* 8 Primeros    vecinos   */
    jp[2] = j + neigh_p[2];
    jm[2] = j + neigh_n[2];
    jp[3] = j + neigh_p[3];
    jm[3] = j + neigh_n[3];

}




void Genera_Conf(int flag)
{
      
       o4v temp,uno,cero,fijo;
   
       int igen,x,y,z,t,mu;
       float invnorma2,norma2;
       float a,b,c,d;
       uno.a1=uno.a2=uno.a3=0.0F;
       uno.a0=1.0F;

       fijo.a0=fijo.a1=0.5;
       fijo.a2=1.0F/(2.0*sqrt(2.0));
       fijo.a3=sqrt(3.0F)/(2.0*sqrt(2.0));
       
        
       cero = uno;
       cero.a0 = 0.0F;
       
       igen=0;
  

       
       if(flag < 2)
       {
           for(t=0;t<L;t++)
               for(z=0;z<L;z++)
                   for(y=0;y<L;y++)
                       for(x=0;x<L;x++)
                       {
                                    
                           if (flag == 0)       /* Config  random  */
                           { 


                               /* determinista */
                               
                               a=FRANDOM*2.0F - 1.0F; 
                               b=FRANDOM*2.0F - 1.0F;
                               c=FRANDOM*2.0F - 1.0F;
                               d=FRANDOM*2.0F - 1.0F;
                               norma2=a*a + b*b + c*c + d*d;
                           
                               
                               temp.a0=a;
                               temp.a1=b;
                               temp.a2=c;
                               temp.a3=d;
                               
                               invnorma2=1.0F/sqrt((double)norma2);
                               
                               _multesc(temp,temp,invnorma2);
                               
                               
                               
                               f[igen]=temp;  /* temp */
                           
                               for(mu = 0; mu<4; mu++)
                               {
                                   
                                   a=FRANDOM*2.0F - 1.0F; 
                                   b=FRANDOM*2.0F - 1.0F;
                                   c=FRANDOM*2.0F - 1.0F;
                                   d=FRANDOM*2.0F - 1.0F;
                                   norma2=a*a + b*b + c*c + d*d;
                                   temp.a0=a;
                                   temp.a1=b;
                                   temp.a2=c;
                                   temp.a3=d;
                                   
                                   invnorma2=1.0F/sqrt((double)norma2);
                                   
                                   _multesc(temp,temp,invnorma2);
                                   
                                   U[mu][igen] = temp;  /* temp */
                               
                               
                               }
                           
                           }                  

                           else
                           {
                               if(flag==1)
                               {   
                                
                                   
                                   f[igen] = uno;          

                                   for(mu=0;mu<4;mu++) U[mu][igen] = uno;
                               }
                           } 
#ifdef PUROGAUGE

f[igen] = uno;
  
#endif                           
                           igen++;
                           
                       }
           
       }
       
       

        
   }

void Lee_Datos(void)
{
    

    int j;
    int * ptdatos_int;
    float * ptdatos_real;
    char name[33];
       

    Finput=fopen("input","r");
    if (Finput==NULL)
    {
        printf(" No existe el fichero 'input'.\n");
        exit(0);
    }
     

    fscanf(Finput,"%s",dir);  

    for (j=0,ptdatos_int=&datos.itmax;j<NDATINT;j++)
        fscanf(Finput,"%d",ptdatos_int++);
    for (j=0,ptdatos_real=&datos.kappa1;j<NDATFLOAT;j++)
        fscanf(Finput,"%f",ptdatos_real++);
    fclose(Finput);
#ifdef DEBUG
pinta_datos(&datos);
#endif
    if(datos.flag < 2)
    {
       sprintf(name,"%s%s%03d.DAT",dir,"OUT",datos.itcut);

        Foutput=fopen(name,"rb");
        if (Foutput!=NULL)
        {
            fclose(Foutput);
            printf(" %s  ya existe.\a\n",name);
            exit;
        }

       sprintf(name,"%s%s",dir,"conf");

        Foutput=fopen(name,"rb");
        if (Foutput!=NULL)
        {
            fclose(Foutput);
            printf(" %s  ya existe.\a\n",name);
            exit;
        }
    }

}



void Lee_Conf(int i)
{
    struct s_datos datosb;
    char name[LPATH+33];
    int j;    
    

    sprintf(name,"%s%s",dir,"conf");

    Fconfig=fopen(name,"rb");
    if (Fconfig==NULL)
    {
        printf(" No existe el fichero '%s'.\a\n",name);
        exit;
    }

    fread(&datosb,sizeof(datosb),1,Fconfig);
 
    fread(f,sizeof(o4v),V,Fconfig);

    fread(U,sizeof(o4v),4*V,Fconfig);
    
#ifdef DEBUG
    pinta_datos(&datosb);
#endif

    if (datos.itmax   != datosb.itmax ||
        datos.mfresh != datosb.mfresh ||
        datos.kappa1   != datosb.kappa1 ||
        datos.kappa2   != datosb.kappa2 ||
        datos.beta     != datosb.beta ||
        datos.flag  == 3 )
    {
        printf("Datos de '%s' son incompatibles con los de 'input'\n",name);
        printf("Se usaran los datos de `input` para simular\n");

        sprintf(name,"%s%s%03d.DAT",dir,"OUT",datos.itcut);

        Foutput=fopen(name,"rb");
        if (Foutput!=NULL)
        {
            fclose(Foutput);
            printf(" %s  ya existe.\a\n",name);
            exit;
        }
    }
    else
    {
        datos.itcut=datosb.itcut;
        datos.seed=datosb.seed;
    }
    fclose(Fconfig);
}


void escribe_conf(int i)
{
    char name[LPATH+33],name_dollar[LPATH+33],name_old[LPATH+33];


    sprintf(name_dollar,"%s%s",dir,"conf.$$$");
    sprintf(name,"%s%s",dir,"conf");
    sprintf(name_old,"%s%s",dir,"conf.old");

    Fconfig=fopen(name_dollar,"wb");
    fwrite(&datos,sizeof(datos),1,Fconfig);
#ifdef DEBUG
    pinta_datos(&datos);
#endif
   
    fwrite(f,sizeof(o4v),V,Fconfig);
    fwrite(U,sizeof(o4v),4*V,Fconfig);
    
    fclose(Fconfig);
    remove(name_old);
    rename(name,name_old);
    rename(name_dollar,name);
}


int pinta_datos(struct s_datos *dat)   
{
    printf("itmax  %d \n",dat->itmax);
    printf("mfresh  %d \n",dat->mfresh);
    printf("nbin   %d \n",dat->nbin);
    printf("itcut  %d \n",dat->itcut);
    printf("flag   %d \n",dat->flag);
    printf("seed   %d \n",dat->seed);
    printf("kappa1  %f \n",dat->kappa1);
    printf("kappa2  %f \n",dat->kappa2);
    printf("beta    %f \n",dat->beta);
#ifdef HISTER
    printf("dbeta %f \n",dat->dbeta);
    printf("dkappa2 %f \n",dat->dkappa2);
#endif
    
    return(0);
}         


void Inicializa(int flag)
{

   if( flag  >= 2) Lee_Conf(0);  

   if( flag < 2 )  Genera_Conf(flag);
     
 
}


void Escribe_Result( int i)
{
    int idat;
    char name[LPATH+33];

    sprintf(name,"%s%s%03d.DAT",dir,"OUT",i);
    
    Foutput=fopen(name,"wb");

    fwrite(&datos,sizeof(datos),1,Foutput);
    for(idat=0;idat<n_obs;idat++)   
       fwrite(&obs_dat[idat][0],4u*datos.itmax,1,Foutput);

    fclose(Foutput);


}


void escribe_hister(int i)
{
#ifdef HISTER
    FILE *Fhister;
    char name[LPATH+33];
    int k;


    sprintf(name,"%s%s.dat",dir,"hister");

    Fhister=fopen(name,"a");    
    fprintf(Fhister,"%7.4f %7.4f ",datos.kappa2,datos.beta);
    for (k=0;k<NOBS_HISTER;k++)
        fprintf(Fhister,"%f ",resultados[k]);
    fprintf(Fhister,"\n");    
    fclose(Fhister);   

    printf("(%5.3f,%5.3f),",datos.kappa2,datos.beta);

#endif /* HISTER */

    printf("i=%3d,E_g=%7.4f,",i,
            resultados[0]);
    printf("E_1=%7.4f,",
            resultados[1]);
    printf("E_2=%7.4f,",
            resultados[2]);


}



