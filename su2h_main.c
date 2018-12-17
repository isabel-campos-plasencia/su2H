

/* MAIN PROGRAM FOR SU2-HIGGS */

#define MAIN
#include "su2h.h"
#include "su2h_arit.h"


/* DECLARACIONES DE FUNCIONES EXTERNAS AL MAIN */

extern void Init_Rand(int); 
extern void Overrelaxed(int); 
extern void Direccionamientos(void);
extern void Lee_Datos(void);
extern void Genera_Conf(int);
extern void Inicializa(int);
extern void Heat_Bath(int); 
extern void StapleLinks(int,int);
extern void StapleHiggs(int);
extern void Energias(void);
extern void Energia_g(int);
extern void Energia_1(int);
extern void Energia_2(int); 
extern void Neigh(int);
extern void Trans_Gauge(void);
extern void Energia_update(void);
extern void escribe_hister(int);
extern void Vecinos(int);
extern void Links_2_vecinos(int);
extern void Order(void);
extern void Magnet(void);
extern void SD(void);


o4v f[V],U[4][V];

o4v stapleH,sG,sg2[4];

o4v fi,fia,sW,sW1[4],sW2[4];

float resultados[NOBS_HISTER];

float kappa1,kappa2,beta;

float e_g,e_1,e_2,E_g,E_1,E_2,SDW,e1;

int x_p[L],y_p[L],z_p[L],t_p[L],
    x_m[L],y_m[L],z_m[L],t_m[L];

int seed,flag,site;

int jpx,jpy,jpz,jpt,jmx,jmy,jmz,jmt;
int jpp[6],jpm[6],jmp[6],jmm[6];
int neigh_px,neigh_py,neigh_pz,neigh_pt;
int neigh_mx,neigh_my,neigh_mz,neigh_mt;

int neigh_p[4],neigh_n[4],jp[4],jm[4];


struct s_datos datos;

float obs_dat[n_obs][maxit];


int main(void)
{
    
    int x,y,z,t,mfresh,it,mfr,ibin,i,nodo,vez;
    float m;
    o4v magnet[16],cero;
    cero.a0 =cero.a1=cero.a2=cero.a3=0.0F;

    Direccionamientos();

    
    Lee_Datos();


    Init_Rand(datos.seed);
    

    Inicializa(datos.flag);

    mfresh = datos.mfresh;     

    for(ibin=datos.itcut;ibin<datos.nbin;ibin++)
    {
        Init_Rand((unsigned) datos.seed);
    


#ifdef HISTER

	

	if(ibin<datos.nbin/2)
	  {
	    datos.beta+=datos.dbeta;
	    datos.kappa2+=datos.dkappa2;
	    
	  }
	else
	  {
	    datos.beta-=datos.dbeta;
	    datos.kappa2-=datos.dkappa2;
	  }
    
#endif
        
	kappa1 = datos.kappa1;
	kappa2 = datos.kappa2;
	beta = datos.beta;
	

        /*  programa principal  */
    
        for(it=0;it<datos.itmax;it++)
        {
            for(mfr=0;mfr<mfresh;mfr++)
            {
                site = 0;
            
                for(t=0;t<L;t++)
                {        
                    neigh_pt = neigh_p[3] = t_p[t];
                    neigh_mt = neigh_n[3] = t_m[t];
                    
                    for(z=0;z<L;z++)
                    {
                       neigh_pz = neigh_p[2] = z_p[z];
                       neigh_mz = neigh_n[2] = z_m[z];
                        
                        for(y=0;y<L;y++)
                        {
                           neigh_py = neigh_p[1] = y_p[y];
                           neigh_my = neigh_n[1] = y_m[y];
                            
                            for(x=0;x<L;x++)
                            {
                               neigh_px = neigh_p[0] = x_p[x];
                               neigh_mx = neigh_n[0] = x_m[x];
                                
			       /* if(mfr == 0) */   Heat_Bath(site);  
			       /* else    Overrelaxed(site); */     
				
                   
                                site++;
                            
                            }   /* coord x  */
                        }       /* coord y  */
                    }           /* coord z  */
                }               /* coord t  */
            }                   /* fin de mfresh */
            

            /* MEDIDAS  */

            Energias();
          
            obs_dat[0][it] = E_g * Normaenerp;
            obs_dat[1][it] = E_1 * Normaener1;
            obs_dat[2][it] = E_2 * Normaener2;


        }      /* fin de it */
    
	

#ifdef HISTER

      resultados[0]=obs_dat[0][datos.itmax-1];
      resultados[1]=obs_dat[1][datos.itmax-1];
      resultados[2]=obs_dat[2][datos.itmax-1];
      
      escribe_hister(ibin);  /* fichero de histeresis */

#endif

#ifndef HISTER    

        Escribe_Result(ibin);

	printf("Energia plaquetas = %f\n",obs_dat[0][datos.itmax-1]);
	printf("Energia 1os vecinos = %f\n",obs_dat[1][datos.itmax-1]);
	printf("Energia 2os vecinos = %f\n",obs_dat[2][datos.itmax-1]);

        
#endif
        datos.seed = RANDOM;
        datos.itcut = ibin + 1;
    
    
        escribe_conf(0);
    

    }   /*******  FINAL DE NBIN  **********/


    return 0;
}
