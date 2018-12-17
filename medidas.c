
/* FUNCIONES DE MEDIDA DE OBSERVABLES  SU(2)-HIGGS  */

#include "su2h.h"
#include "su2h_arit.h"


extern o4v f[V],U[4][V],stapleH,sg2[4],sW,sW1[],sW2[];

extern int neigh_p[4],neigh_n[4],x_p[],y_p[],z_p[],t_p[],jm[];
extern int jp[],jm[],x_m[],y_m[],z_m[],t_m[];

extern int jpx,jpy,jpz,jpt,jmx,jmy,jmz,jmt;
extern int jpp[],jpm[],jmp[],jmm[];
extern int neigh_px,neigh_py,neigh_pz,neigh_pt;
extern int neigh_mx,neigh_my,neigh_mz,neigh_mt;

extern float kappa1,kappa2,beta,e_g,e_1,e_2,E_g,E_1,E_2,SDW,e1;


o4v temp1,temp2,temp3,inter1,inter2,inter3;


void Energia_g(int s)
{
    int mu,nu;

    
    e_g = 0.0;
    
    for(nu=0;nu<4;nu++)
    {
        
        for(mu=0;mu<nu;mu++)
        {
            _multlink(inter1,U[mu][s],U[nu][s+neigh_p[mu]]);
            _adjunta(temp1,U[mu][s+neigh_p[nu]]);
            _adjunta(temp2,U[nu][s]);
            
            _multlink(inter2,temp1,temp2);
            
            _multlink(temp1,inter1,inter2);
            
	    e_g += (1- temp1.a0);
            
            

        }
    }
    
        
}


void Energia_1(int s)
{

    int mu;
    float prueba_e;
    o4v t1;
    
    
    e_1 = 0.0;
    prueba_e = 0.0;
    
    for(mu=0;mu<4;mu++)
    {
        _adjunta(t1,f[s + neigh_p[mu]]);
        _multlink(temp2,f[s],U[mu][s]);
        
        _multlink(temp1,temp2,t1);
    

        e_1 += (temp1.a0);
        

    }
    
        
}  

void Energia_2(int s)
{
    int mu,nu;
    
    e_2 = 0.0;
    
    for(nu=0;nu<4;nu++)
        for(mu=0;mu<nu;mu++)
        {
            _multlink(temp1,U[mu][s],U[nu][s + neigh_p[mu]]);
            _multlink(temp2,U[nu][s],U[mu][s + neigh_p[nu]]);
           
            _adjunta(inter1,f[s]);
            _multlink(inter2,f[s + neigh_p[mu] + neigh_p[nu]],inter1);
            
            _suma(inter1,temp1,temp2);
            
            _multlink(temp1,inter2,inter1);

            
            e_2 += temp1.a0;
        }
    
    e_2 = 0.5 * e_2;

    
}                   


void Energias(void)
{
    int x,y,z,t,s,m;
    float e2,p,cociente;
    
    p = E_g = E_1 =E_2 = 0.0F;
    
    s = 0;
            
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
			    Energia_g(s); 

                            Energia_1(s);   
			    Energia_2(s);   

                            
			    E_g += e_g;  
                            E_1 += e_1; 
			    E_2 += e_2; 
		
                            s++;
                            
                        }   /* coord x  */
                    }       /* coord y  */
                }           /* coord z  */
            }               /* coord t  */


/*    printf("E plaquetas = %f\n",E_g/(6*V));
    printf("E 1os vecinos = %f\n",E_1/(4*V));
    printf("E 2os vecinos = %f\n",E_2/(6*V));
  */ 
}

void SD(void)
{

  o4v cero,W,t1,t2;
  int x,y,z,t,s;
  int mu,nu,m,n;
  float sd,primer,p,e2,coupling,num,den;
    

  s = 0;  
  cero.a0 = cero.a1 = cero.a2 = cero.a3 = 0.0F;
  E_2 = E_1 = 0.0;
  p =  primer = e_2 = sd=0.0;
  num = den = 0.0;

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
		  
		  Vecinos(s);
		  
		  Links_2_vecinos(s);

		  W =cero;

		  for(m=0;m<4;m++)
		    {
		      _adjunta(t1,U[m][s]);
		      
		      W.a0 += t1.a0 + U[m][s + neigh_n[m]].a0;
		      W.a1 += t1.a1 + U[m][s + neigh_n[m]].a1;
		      W.a2 += t1.a2 + U[m][s + neigh_n[m]].a2;
		      W.a3 += t1.a3 + U[m][s + neigh_n[m]].a3;

		    }


		  num += W.a0;
		  den += (W.a1*W.a1 + W.a2*W.a2 + W.a3*W.a3)/3.0;
	      

		  for(m=0;m<4;m++)
		    {
		      StapleLinks(s,m);
		      sd += sW.a1*sW.a1 + sW.a2*sW.a2 + sW.a3*sW.a3;

		      primer += 2.0*kappa1*
			(U[m][s].a1*sW.a1 + U[m][s].a2*sW.a2 
			 + U[m][s].a3*sW.a3);
		      e2 += (sW.a0);
		    }

	
		  s++;
		  
		}
	    }
	}
    }

  printf(" beta = %f\n",((3.0*e2) - primer)/sd);  

  printf("kappa = %f\n",num/(2.0*den));  
    

}



