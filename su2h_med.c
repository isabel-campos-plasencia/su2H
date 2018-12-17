
/* FUNCIONES DE MEDIDA DE OBSERVABLES  SU(2)-HIGGS  */

#include "su2h.h"
#include "su2h_arit.h"


extern o4v stapleH,sG,sW,sW1[4],sW2;

extern o4v f[V],U[4][V];

extern int neigh_p[4],neigh_n[4],x_p[],y_p[],z_p[],t_p[],jm[];
extern int jp[],jm[],x_m[],y_m[],z_m[],t_m[];

extern int jpx,jpy,jpz,jpt,jmx,jmy,jmz,jmt;
extern int jpp[],jpm[],jmp[],jmm[];
extern int neigh_px,neigh_py,neigh_pz,neigh_pt;
extern int neigh_mx,neigh_my,neigh_mz,neigh_mt;

extern float kappa1,kappa2,beta,e_g,e_1,e_2,E_g,E_1,E_2,SDW,e1;


o4v nn[24],temp1,temp2,temp3,inter1,inter2,inter3;


void Vecinos(int j)
{

  jp[0] = jpx = j + neigh_px;        /* CALCULO DE LOS SITES VECINOS */
  jm[0] = jmx = j + neigh_mx;
  jp[1] = jpy = j + neigh_py;
  jm[1] = jmy = j + neigh_my;        /* 8 Primeros    vecinos   */
  jp[2] = jpz = j + neigh_pz;
  jm[2] = jmz = j + neigh_mz;
  jp[3] = jpt = j + neigh_pt;
  jm[3] = jmt = j + neigh_mt;

                             /* 24 segundos vecinos */
  jpp[0] = jpx + neigh_py;
  jpp[1] = jpx + neigh_pz;
  jpp[2] = jpx + neigh_pt;
  jpp[3] = jpy + neigh_pz;
  jpp[4] = jpy + neigh_pt;
  jpp[5] = jpz + neigh_pt;

  jpm[0] = jpx + neigh_my;
  jpm[1] = jpx + neigh_mz;
  jpm[2] = jpx + neigh_mt;
  jpm[3] = jpy + neigh_mz;
  jpm[4] = jpy + neigh_mt;
  jpm[5] = jpz + neigh_mt;

  jmp[0] = jmx + neigh_py;
  jmp[1] = jmx + neigh_pz;
  jmp[2] = jmx + neigh_pt;
  jmp[3] = jmy + neigh_pz;
  jmp[4] = jmy + neigh_pt;
  jmp[5] = jmz + neigh_pt;

  jmm[0] = jmx + neigh_my;
  jmm[1] = jmx + neigh_mz;
  jmm[2] = jmx + neigh_mt;
  jmm[3] = jmy + neigh_mz;
  jmm[4] = jmy + neigh_mt;
  jmm[5] = jmz + neigh_mt;


}






void StapleHiggs(int s)
{
    int m,n,i;
    o4v t1,t2,t3,t4,t5,t6,inter1,inter2,inter3,inter4;
    o4v cero,s2,s1,fh,fha;
    
    
    cero.a0 = cero.a1 = cero.a2 = cero.a3 = 0.0F;
    
    fh = f[s];
    _adjunta(fha,fh);

    s1 = cero;
    s2 = cero;

    Neigh(s);
    

    for(m=0;m<4;m++)           
    {
        
        _multlink_t(t1,U[m][s],f[jp[m]]); 
        
	_multlink_t_t(t3,U[m][jm[m]],f[jm[m]]); 
	
        
        s1.a0 += (t1.a0 + t3.a0);
        s1.a1 += (t1.a1 + t3.a1);
        s1.a2 += (t1.a2 + t3.a2);
        s1.a3 += (t1.a3 + t3.a3);
     

    }
    

    _multesc(s1,s1,2.0F*kappa1);  /* LA TRAZA ES 2 * s1.a0 */
     
    /*** SEGUNDOS  VECINOS + Links  ***/
 
    /* jpp segundos ++ */

    Vecinos(s);


    _multlink(t1,U[1][s],U[0][jpy]);
    _multlink(t2,U[0][s],U[1][jpx]);
    _suma(t3,t1,t2);
    _multlink(nn[0],t3,f[jpp[0]]);


    _multlink(t1,U[2][s],U[0][jpz]);
    _multlink(t2,U[0][s],U[2][jpx]);
    _suma(t3,t1,t2);
    _multlink(nn[1],t3,f[jpp[1]]);


    _multlink(t1,U[3][s],U[0][jpt]);
    _multlink(t2,U[0][s],U[3][jpx]);
    _suma(t3,t1,t2);
    _multlink(nn[2],t3,f[jpp[2]]);


    _multlink(t1,U[2][s],U[1][jpz]);
    _multlink(t2,U[1][s],U[2][jpy]);
    _suma(t3,t1,t2);
    _multlink(nn[3],t3,f[jpp[3]]);


    _multlink(t1,U[3][s],U[1][jpt]);
    _multlink(t2,U[1][s],U[3][jpy]);
    _suma(t3,t1,t2);
    _multlink(nn[4],t3,f[jpp[4]]);


    _multlink(t1,U[3][s],U[2][jpt]);
    _multlink(t2,U[2][s],U[3][jpz]);
    _suma(t3,t1,t2);
    _multlink(nn[5],t3,f[jpp[5]]);

    /* jpm */

    _adjunta(t1,U[1][jpm[0]]);
    _adjunta(t2,U[1][jmy]);
    _multlink(t3,U[0][s],t1);
    _multlink(t4,t2,U[0][jmy]);
    _suma(t5,t3,t4);
    _multlink(nn[6],t5,f[jpm[0]]);


    _adjunta(t1,U[2][jpm[1]]);
    _adjunta(t2,U[2][jmz]);
    _multlink(t3,U[0][s],t1);
    _multlink(t4,t2,U[0][jmz]);
    _suma(t5,t3,t4);
    _multlink(nn[7],t5,f[jpm[1]]);


    _adjunta(t1,U[3][jpm[2]]);
    _adjunta(t2,U[3][jmt]);
    _multlink(t3,U[0][s],t1);
    _multlink(t4,t2,U[0][jmt]);
    _suma(t5,t3,t4);
    _multlink(nn[8],t5,f[jpm[2]]);


    _adjunta(t1,U[2][jpm[3]]);
    _adjunta(t2,U[2][jmz]);
    _multlink(t3,U[1][s],t1);
    _multlink(t4,t2,U[1][jmz]);
    _suma(t5,t3,t4);
    _multlink(nn[9],t5,f[jpm[3]]);    


    _adjunta(t1,U[3][jpm[4]]);
    _adjunta(t2,U[3][jmt]);
    _multlink(t3,U[1][s],t1);
    _multlink(t4,t2,U[1][jmt]);
    _suma(t5,t3,t4);
    _multlink(nn[10],t5,f[jpm[4]]);    


    _adjunta(t1,U[3][jpm[5]]);
    _adjunta(t2,U[3][jmt]);
    _multlink(t3,U[2][s],t1);
    _multlink(t4,t2,U[2][jmt]);
    _suma(t5,t3,t4);
    _multlink(nn[11],t5,f[jpm[5]]);    


    /* jmp */

    _adjunta(t1,U[0][jmp[0]]);
    _adjunta(t2,U[0][jmx]);
    _multlink(t3,U[1][s],t1);
    _multlink(t4,t2,U[1][jmx]);
    _suma(t5,t3,t4);
    _multlink(nn[12],t5,f[jmp[0]]);

    _adjunta(t1,U[0][jmp[1]]);
    _adjunta(t2,U[0][jmx]);
    _multlink(t3,U[2][s],t1);
    _multlink(t4,t2,U[2][jmx]);
    _suma(t5,t3,t4);
    _multlink(nn[13],t5,f[jmp[1]]);

    _adjunta(t1,U[0][jmp[2]]);
    _adjunta(t2,U[0][jmx]);
    _multlink(t3,U[3][s],t1);
    _multlink(t4,t2,U[3][jmx]);
    _suma(t5,t3,t4);
    _multlink(nn[14],t5,f[jmp[2]]);

    _adjunta(t1,U[1][jmp[3]]);
    _adjunta(t2,U[1][jmy]);
    _multlink(t3,U[2][s],t1);
    _multlink(t4,t2,U[2][jmy]);
    _suma(t5,t3,t4);
    _multlink(nn[15],t5,f[jmp[3]]);

    _adjunta(t1,U[1][jmp[4]]);
    _adjunta(t2,U[1][jmy]);
    _multlink(t3,U[3][s],t1);
    _multlink(t4,t2,U[3][jmy]);
    _suma(t5,t3,t4);
    _multlink(nn[16],t5,f[jmp[4]]);

    _adjunta(t1,U[2][jmp[5]]);
    _adjunta(t2,U[2][jmz]);
    _multlink(t3,U[3][s],t1);
    _multlink(t4,t2,U[3][jmz]);
    _suma(t5,t3,t4);
    _multlink(nn[17],t5,f[jmp[5]]);


    /* jmm */

    _multlink(t1,U[1][jmm[0]],U[0][jmx]);
    _multlink(t2,U[0][jmm[0]],U[1][jmy]);
    _suma(t3,t1,t2);
    _adjunta(t3,t3);
    _multlink(nn[18],t3,f[jmm[0]]);

    _multlink(t1,U[2][jmm[1]],U[0][jmx]);
    _multlink(t2,U[0][jmm[1]],U[2][jmz]);
    _suma(t3,t1,t2);
    _adjunta(t3,t3);
    _multlink(nn[19],t3,f[jmm[1]]);

    _multlink(t1,U[3][jmm[2]],U[0][jmx]);
    _multlink(t2,U[0][jmm[2]],U[3][jmt]);
    _suma(t3,t1,t2);
    _adjunta(t3,t3);
    _multlink(nn[20],t3,f[jmm[2]]);

    _multlink(t1,U[2][jmm[3]],U[1][jmy]);
    _multlink(t2,U[1][jmm[3]],U[2][jmz]);
    _suma(t4,t1,t2);
    _adjunta(t3,t4);
    _multlink(nn[21],t3,f[jmm[3]]);

    _multlink(t1,U[3][jmm[4]],U[1][jmy]);
    _multlink(t2,U[1][jmm[4]],U[3][jmt]);
    _suma(t3,t1,t2);
    _adjunta(t3,t3);
    _multlink(nn[22],t3,f[jmm[4]]);

    _multlink(t1,U[3][jmm[5]],U[2][jmz]);
    _multlink(t2,U[2][jmm[5]],U[3][jmt]);
    _suma(t3,t1,t2);
    _adjunta(t3,t3);
    _multlink(nn[23],t3,f[jmm[5]]);
    
    
    /*******/
    
    
    for(i=0;i<24;i++)
      {
	s2.a0 += nn[i].a0;
	s2.a1 += nn[i].a1;
	s2.a2 += nn[i].a2;
	s2.a3 += nn[i].a3;
      }

    _adjunta(s2,s2);
 

    _multesc(s2,s2,kappa2);

 
    _suma(stapleH,s1,s2); 


}




void StapleLinks(int s, int m) 
{
    int n;
    o4v cero,sgg,sg2,sg1,uno,A;
    float den,A0;
    o4v t1,t2,t3,t4,t5,fia;
    
    
    cero.a0 = cero.a1 = cero.a2 = cero.a3 = 0.0F;
    sgg = cero;
    sg2 = cero;
    sg1 = cero;
    
    uno = cero;
    uno.a0 = 1.0F;
    


#ifndef PUROGAUGE

    _adjunta(fia,f[jp[m]]);
    
    _multlink(sg1,fia,f[s]);
    
#endif
    
    for(n=0;n<4;n++)
    {
        if(n != m)
        {
            _adjunta(temp1,U[m][jp[n]]);
            _multlink(temp2,U[n][jp[m]],temp1);
            _adjunta(temp3,U[n][s]);
            
            _multlink(inter1,temp2,temp3);
            
            _adjunta(temp1,U[n][jp[m] + neigh_n[n]]);
            _adjunta(temp2,U[m][jm[n]]);
            _multlink(temp3,temp1,temp2);
            _multlink(inter2,temp3,U[n][jm[n]]);
 
            sgg.a0 += (inter1.a0 + inter2.a0);
            sgg.a1 += (inter1.a1 + inter2.a1);
            sgg.a2 += (inter1.a2 + inter2.a2);
            sgg.a3 += (inter1.a3 + inter2.a3);
            
        }
    }


    _multesc(sg1,sg1,2.0F*kappa1);  /* LA TRAZA ES 2*sg1.a0 */

    _multlink(sW,U[m][s],sgg); /* variable de S-D */

    _multesc(sgg,sgg,beta);   
    
    _suma(sG,sgg,sg1); 
    

}





void Trans_Gauge(void)
{
    int x1,x2,x3,x4,igen,mu;
    float a,b,c,d,norma2,invnorma2;
    
    
    o4v temp,inter,g[V];

    igen=0;
    
    for(x4=0;x4<L;x4++)
               for(x3=0;x3<L;x3++)
                   for(x2=0;x2<L;x2++)
                       for(x1=0;x1<L;x1++)
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
                           
                           g[igen] = temp;

                           igen++;
                           
                       }
    


    igen = 0;
    
    Energias();
    

    for(x4=0;x4<L;x4++)
    {
        neigh_p[3] = t_p[x4];
        neigh_n[3] = t_m[x4];
        
        for(x3=0;x3<L;x3++)
        {
            neigh_p[2] = z_p[x3];
            neigh_n[2] = z_m[x3];
            
            for(x2=0;x2<L;x2++)
            {
                neigh_p[1] = y_p[x2];
                neigh_n[1] = y_m[x2];
                
                for(x1=0;x1<L;x1++)
                {
                    neigh_p[0] = x_p[x1];
                    neigh_n[0] = x_m[x1];
                    
                    Neigh(igen);
                    
                    for(mu=0;mu<4;mu++)
                    {
                        _adjunta(temp,g[jp[mu]]);
                        _multlink(inter,U[mu][igen],temp);
                        _multlink(U[mu][igen],g[igen],inter);
                    }
                    igen++;
                    
                }
            }
        }
    }
    
    printf("Transformacion Gauge \n");

    Energias();
    



}


void Energia_update(void)
{
    
   int x,y,z,t,s,n,m;
   o4v temp1,temp2,temp3,inter1,inter2,inter3,sgg,cero;
   
   cero.a0=cero.a1=cero.a2=cero.a3=0.0F;
   
    E_g = E_1 =E_2 = 0.0F;
    
    s = 0;
            
            for(t=0;t<L;t++)
            {        
                neigh_p[3] = t_p[t];
                neigh_n[3] = t_m[t];
                
                for(z=0;z<L;z++)
                {
                    neigh_p[2] = z_p[z];
                    neigh_n[2] = z_m[z];
                    
                    for(y=0;y<L;y++)
                    {
                        neigh_p[1] = y_p[y];
                        neigh_n[1] = y_m[y];
                        
                        for(x=0;x<L;x++)
                        {
                            neigh_p[0] = x_p[x];
                            neigh_n[0] = x_m[x];
                            
                            for(m=0;m<4;m++)
                            {

                                sgg = cero;
                                Neigh(s);
 
                                for(n=0;n<4;n++)
                                {
                                    if(n != m)
                                    {
                                        _adjunta(temp1,U[m][jp[n]]);
                                        _multlink(temp2,U[n][jp[m]],temp1);
                                        _adjunta(temp3,U[n][s]);
                                        
                                        _multlink(inter1,temp2,temp3);
                        
                                        
                                        
                                        _adjunta(temp1,U[n][jp[m] + neigh_n[n]]);
                                        _adjunta(temp2,U[m][jm[n]]);
                                        _multlink(temp3,temp1,temp2);
                                        _multlink(inter2,temp3,U[n][jm[n]]);
                                        
                                        sgg.a0 += (inter1.a0 + inter2.a0);
                                        sgg.a1 += (inter1.a1 + inter2.a1);

                                        sgg.a2 += (inter1.a2 + inter2.a2);
                                        sgg.a3 += (inter1.a3 + inter2.a3);
                            
                                        
                                    }
                                    
                                    
                                }


                                
                                _multlink(temp1,sgg,U[m][s]);
                                
                                E_g +=(temp1.a0)/4.0F;
                                

                            }
                            
                            s++;
                            
                        }   /* coord x  */
                    }       /* coord y  */
                }           /* coord z  */
            }               /* coord t  */


    printf("E plaquetas update = %f\n",E_g/(6*V));
 
}





    
