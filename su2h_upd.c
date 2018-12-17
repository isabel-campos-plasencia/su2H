
/*  UPDATE  PARA  SU(2)-HIGGS   */

#include "su2h.h"
#include "su2h_arit.h"

extern o4v f[V],U[4][V],stapleH,sG,sg2[4];

extern float kappa1,kappa2,beta;

extern int neigh_p[],neigh_n[];

extern int jpx,jpy,jpz,jpt,jmx,jmy,jmz,jmt;
extern int jpp[],jpm[],jmp[],jmm[];
extern int neigh_px,neigh_py,neigh_pz,neigh_pt;
extern int neigh_mx,neigh_my,neigh_mz,neigh_mt;



Links_2_vecinos(int s)
{
  o4v ss[7],t1,t2,t3,cero,fia;
  int i;
  o4v v1,v2,v3,v4,v5,v6;
  o4v p1,p2,p3,p4,p5,p6;

  cero.a0 = cero.a1 = cero.a2 =cero.a3 =0.0F;
  t1=cero;

  _adjunta(t1,U[1][jpx + neigh_my]);         /* LINK 0 */
  _adjunta(t2,U[2][jpx + neigh_mz]);
  _adjunta(t3,U[3][jpx + neigh_mt]);

  _multlink(ss[0],U[1][jpx],f[jpx + neigh_py]);
  _multlink(ss[1],t1,f[jpx + neigh_my]);
  _multlink(ss[2],U[2][jpx],f[jpx + neigh_pz]);
  _multlink(ss[3],t2,f[jpx + neigh_mz]);
  _multlink(ss[4],U[3][jpx],f[jpx + neigh_pt]);
  _multlink(ss[5],t1,f[jpx + neigh_mt]);


  _adjunta(v1,f[jpy]);
  _adjunta(v2,f[jmy]);
  _adjunta(v3,f[jpz]);
  _adjunta(v4,f[jmz]);
  _adjunta(v5,f[jpt]);
  _adjunta(v6,f[jmt]);
  _adjunta(t1,U[1][s]);
  _adjunta(t2,U[2][s]);
  _adjunta(t3,U[3][s]);

  _multlink(p1,v1,t1);
  _multlink(p2,v2,U[1][jmy]);
  _multlink(p3,v3,t2);
  _multlink(p4,v4,U[2][jmz]);
  _multlink(p5,v5,t3);
  _multlink(p6,v6,U[3][jmt]);

  _suma6(v1,p1,p2,p3,p4,p5,p6);

  _multlink(ss[6],f[jpx],v1);



  for(i=0;i<6;i++)
    {
      t1.a0 += ss[i].a0;
      t1.a1 += ss[i].a1;
      t1.a2 += ss[i].a2;
      t1.a3 += ss[i].a3;
    }

  _adjunta(fia,f[s]);
  _multlink(t2,t1,fia);
  _suma(sg2[0],t2,ss[6]);
  _multesc(sg2[0],sg2[0],kappa2);


  t1 = cero;

  _adjunta(t1,U[0][jpy + neigh_mx]);         /* LINK 1 */
  _adjunta(t2,U[2][jpy + neigh_mz]);
  _adjunta(t3,U[3][jpy + neigh_mt]);

  _multlink(ss[0],U[0][jpy],f[jpy + neigh_px]);
  _multlink(ss[1],t1,f[jpy + neigh_mx]);
  _multlink(ss[2],U[2][jpy],f[jpy + neigh_pz]);
  _multlink(ss[3],t2,f[jpy + neigh_mz]);
  _multlink(ss[4],U[3][jpy],f[jpy + neigh_pt]);
  _multlink(ss[5],t1,f[jpy + neigh_mt]);

  _adjunta(v1,f[jpx]);
  _adjunta(v2,f[jmx]);
  _adjunta(v3,f[jpz]);
  _adjunta(v4,f[jmz]);
  _adjunta(v5,f[jpt]);
  _adjunta(v6,f[jmt]);
  _adjunta(t1,U[0][s]);
  _adjunta(t2,U[2][s]);
  _adjunta(t3,U[3][s]);

  _multlink(p1,v1,t1);
  _multlink(p2,v2,U[0][jmx]);
  _multlink(p3,v3,t2);
  _multlink(p4,v4,U[2][jmz]);
  _multlink(p5,v5,t3);
  _multlink(p6,v6,U[3][jmt]);

  _suma6(v1,p1,p2,p3,p4,p5,p6);

  _multlink(ss[6],f[jpy],v1);




  for(i=0;i<6;i++)
    {
      t1.a0 += ss[i].a0;
      t1.a1 += ss[i].a1;
      t1.a2 += ss[i].a2;
      t1.a3 += ss[i].a3;
    }
  

  _multlink(t2,t1,fia);
  _suma(sg2[1],t2,ss[6]);
  _multesc(sg2[1],sg2[1],kappa2);


  t1 = cero;

  _adjunta(t1,U[0][jpz + neigh_mx]);         /* LINK 2 */
  _adjunta(t2,U[1][jpz + neigh_my]);
  _adjunta(t3,U[3][jpz + neigh_mt]);

  _multlink(ss[0],U[0][jpz],f[jpz + neigh_px]);
  _multlink(ss[1],t1,f[jpz + neigh_mx]);
  _multlink(ss[2],U[1][jpz],f[jpz + neigh_py]);
  _multlink(ss[3],t2,f[jpz + neigh_my]);
  _multlink(ss[4],U[3][jpz],f[jpz + neigh_pt]);
  _multlink(ss[5],t1,f[jpz + neigh_mt]);

  _adjunta(v1,f[jpy]);
  _adjunta(v2,f[jmy]);
  _adjunta(v3,f[jpx]);
  _adjunta(v4,f[jmx]);
  _adjunta(v5,f[jpt]);
  _adjunta(v6,f[jmt]);
  _adjunta(t1,U[1][s]);
  _adjunta(t2,U[0][s]);
  _adjunta(t3,U[3][s]);

  _multlink(p1,v1,t1);
  _multlink(p2,v2,U[1][jmy]);
  _multlink(p3,v3,t2);
  _multlink(p4,v4,U[0][jmx]);
  _multlink(p5,v5,t3);
  _multlink(p6,v6,U[3][jmt]);

  _suma6(v1,p1,p2,p3,p4,p5,p6);

  _multlink(ss[6],f[jpz],v1);

  for(i=0;i<6;i++)
    {
      t1.a0 += ss[i].a0;
      t1.a1 += ss[i].a1;
      t1.a2 += ss[i].a2;
      t1.a3 += ss[i].a3;
    }
  
      _multlink(t2,t1,fia);
      _suma(sg2[2],t2,ss[6]);
      _multesc(sg2[2],sg2[2],kappa2);


  t1 = cero;

  _adjunta(t1,U[0][jpt + neigh_mx]);         /* LINK 3 */
  _adjunta(t2,U[1][jpt + neigh_my]);
  _adjunta(t3,U[2][jpt + neigh_mz]);

  _multlink(ss[0],U[1][jpt],f[jpt + neigh_px]);
  _multlink(ss[1],t1,f[jpt + neigh_mx]);
  _multlink(ss[2],U[2][jpt],f[jpt + neigh_py]);
  _multlink(ss[3],t2,f[jpt + neigh_my]);
  _multlink(ss[4],U[3][jpt],f[jpt + neigh_pz]);
  _multlink(ss[5],t1,f[jpt + neigh_mz]);


  _adjunta(v1,f[jpy]);
  _adjunta(v2,f[jmy]);
  _adjunta(v3,f[jpz]);
  _adjunta(v4,f[jmz]);
  _adjunta(v5,f[jpx]);
  _adjunta(v6,f[jmx]);
  _adjunta(t1,U[1][s]);
  _adjunta(t2,U[2][s]);
  _adjunta(t3,U[0][s]);

  _multlink(p1,v1,t1);
  _multlink(p2,v2,U[1][jmy]);
  _multlink(p3,v3,t2);
  _multlink(p4,v4,U[2][jmz]);
  _multlink(p5,v5,t3);
  _multlink(p6,v6,U[0][jmx]);

  _suma6(v1,p1,p2,p3,p4,p5,p6);

  _multlink(ss[6],f[jpt],v1);




  for(i=0;i<6;i++)
    {
      t1.a0 += ss[i].a0;
      t1.a1 += ss[i].a1;
      t1.a2 += ss[i].a2;
      t1.a3 += ss[i].a3;
    }


  _multlink(t2,t1,fia);
  _suma(sg2[3],t2,ss[6]);
  _multesc(sg2[3],sg2[3],kappa2);




}




void Heat_Bath(int s)
{
    float umag,umagin,damdum,umaginb,norma2,invnorma;
    float a_0,a_1,a_2,a_3,rad,red,x;
    o4v uint, anew,temp,t,stapleG;
    int m; 
    
    Vecinos(s);

    Links_2_vecinos(s); 


#ifndef PUROGAUGE

    StapleHiggs(s);      


    _adjunta(uint,stapleH);  

    
    umag=sqrt(_norm2(uint)); 

    umagin=1.0F/umag;
    
    do
    {
        if ((x=FRANDOM))   
            a_0=1.0F+log(x)*umagin;
        else
            a_0=1.0F-100*umagin;        
 
        rad=1.0F-a_0*a_0;
       
        x=FRANDOM;
 
        
    }
    while (x*x > rad);
    
    do
    {
        a_1=FRANDOM-0.5F;
        a_2=FRANDOM-0.5F;
        red=a_1*a_1+a_2*a_2;
    }
    while (red > 0.25F);

    x=FRANDOM;
    a_3=sqrt(rad)*(2.0F*x-1.0F);
    damdum=rad/red*4.0F*x*(1.0F-x); 
    umaginb=umagin*sqrt(damdum);

    anew.a0 = a_0*umagin;
    anew.a1 = a_3*umagin;
    anew.a2 = a_2*umaginb;
    anew.a3 = a_1*umaginb;


    _rotacion(temp,anew,uint); 

    f[s] = temp;  
    
#endif


            /*  HEAT BATH PARA LAS LINKS  */


    for(m=0;m<4;m++)
    {

       StapleLinks(s,m);             
   

       _suma(stapleG,sG,sg2[m]);

       _adjunta(uint,stapleG);  
    
      umag=sqrt(_norm2(uint)); 

      umagin=1.0F/umag;
      do
        {
	  if ((x=FRANDOM))        
	    a_0=1.0F+log(x)*umagin;
	  else
	    a_0=1.0F-100*umagin;        
	  
	  rad=1.0F-a_0*a_0;
	  
	  x=FRANDOM;
        }
      while (x*x > rad);
      
        do
	  {
            a_1=FRANDOM-0.5F;
            a_2=FRANDOM-0.5F;
            red=a_1*a_1+a_2*a_2;
	  }
        while (red > 0.25F);
        
        x=FRANDOM;
        a_3=sqrt(rad)*(2.0F*x-1.0F);
        damdum=rad/red*4.0F*x*(1.0F-x); 
        umaginb=umagin*sqrt(damdum);
        
        anew.a0 = a_0*umagin;
        anew.a3 = a_3*umagin;
        anew.a2 = a_2*umaginb;
        anew.a1 = a_1*umaginb;
        
    
        _multlink(temp,anew,uint); 
        
        U[m][s] = temp;  
 
    }  
    
 

}  
    





void Overrelaxed(int j)
{
    o4v fi_over,temp,proy,dif,dif_2,U_over,t1;

    o4v stap,uint,anew,stapleG;
    
    float factor,norm_temp,umagin,umag; 
    
    int m;
    
    Vecinos(j);
    Links_2_vecinos(j);
    
#ifndef PUROGAUGE
   
    StapleHiggs(j);

    temp = f[j];
    
    _adjunta(stap,stapleH);  
    
    
    factor = 2.0F*_prodescalar(stap,temp)/_norm2(stap);
    
    fi_over.a0 = stap.a0 * factor - temp.a0;
    fi_over.a1 = stap.a1 * factor - temp.a1;
    fi_over.a2 = stap.a2 * factor - temp.a2;
    fi_over.a3 = stap.a3 * factor - temp.a3;
        
        
    f[j] = fi_over;    
    
    

#endif

    for(m=0;m<4;m++)
    {
        StapleLinks(j,m);
        
        temp = U[m][j];

	_suma(stapleG,sg2[m],sG);

	_adjunta(stap,stapleG);  
    

        factor = 2.0F*_prodescalar(stap,temp)/_norm2(stap);
        
        U_over.a0 = stap.a0 * factor - temp.a0;
        U_over.a1 = stap.a1 * factor - temp.a1;
        U_over.a2 = stap.a2 * factor - temp.a2;
        U_over.a3 = stap.a3 * factor - temp.a3;
        
        U[m][j] = U_over;
        


    
    }
    

}



