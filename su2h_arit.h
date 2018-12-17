

/*  DEFINICION DE LAS OPERACIONES MAS HABITUALES */

#define _multlink(w,v,u) \
	w.a0 = v.a0 * u.a0 - v.a1 * u.a1 - v.a2 * u.a2 - v.a3 * u.a3; \
	w.a1 = v.a0 * u.a1 + v.a3 * u.a2 + v.a1 * u.a0 - v.a2 * u.a3; \
	w.a2 = v.a0 * u.a2 - v.a3 * u.a1 + v.a2 * u.a0 + v.a1 * u.a3; \
	w.a3 = v.a0 * u.a3 + v.a3 * u.a0 + v.a2 * u.a1 - v.a1 * u.a2; 


#define _adjunta(w,u) \
w.a0 =  u.a0; \
w.a1 = -u.a1; \
w.a2 = -u.a2; \
w.a3 = -u.a3;

#define _multlink_t(w,v,u) \
w.a0 = v.a0*u.a0 + v.a1*u.a1 + v.a2*u.a2 + v.a3*u.a3; \
w.a1 = v.a1*u.a0 - v.a0*u.a1 + v.a2*u.a3 - v.a3*u.a2; \
w.a2 = v.a2*u.a0 - v.a0*u.a2 + v.a3*u.a1 - v.a1*u.a3; \
w.a3 = v.a3*u.a0 - v.a0*u.a3 + v.a1*u.a2 - v.a2*u.a1; 

#define _multlink_t_t(w,v,u) \
_multlink(tirar,u,v); \
w.a0 =  tirar.a0; \
w.a1 = -tirar.a1; \
w.a2 = -tirar.a2; \
w.a3 = -tirar.a3;



#define _suma(s,u,v) \
              s.a0 = u.a0 + v.a0; \
              s.a1 = u.a1 + v.a1; \
              s.a2 = u.a2 + v.a2; \
              s.a3 = u.a3 + v.a3;

#define _resta(r,u,v) \
              r.a0 = u.a0 - v.a0; \
              r.a1 = u.a1 - v.a1; \
              r.a2 = u.a2 - v.a2; \
              r.a3 = u.a3 - v.a3;



#define _suma6(s,v1,v2,v3,v4,v5,v6) \
         s.a0 = v1.a0 + v2.a0 + v3.a0 + v4.a0 + v5.a0 + v6.a0; \
         s.a1 = v1.a1 + v2.a1 + v3.a1 + v4.a1 + v5.a1 + v6.a1; \
         s.a2 = v1.a2 + v2.a2 + v3.a2 + v4.a2 + v5.a2 + v6.a2; \
         s.a3 = v1.a3 + v2.a3 + v3.a3 + v4.a3 + v5.a3 + v6.a3; 


#define _suma4(s,v1,v2,v3,v4) \
         s.a0 = v1.a0 + v2.a0 + v3.a0 + v4.a0; \
         s.a1 = v1.a1 + v2.a1 + v3.a1 + v4.a1; \
         s.a2 = v1.a2 + v2.a2 + v3.a2 + v4.a2; \
         s.a3 = v1.a3 + v2.a3 + v3.a3 + v4.a3; 


#define _norm2(u) (u.a0*u.a0 + u.a1*u.a1 + u.a2*u.a2 + u.a3*u.a3)

#define _multesc(u,v,a) \
                  u.a0 = v.a0 *a; \
                  u.a1 = v.a1 *a; \
                  u.a2 = v.a2 *a; \
                  u.a3 = v.a3 *a; 

#define _prodescalar(u,v) (u.a0*v.a0 + u.a1*v.a1 + u.a2*v.a2 + u.a3*v.a3)


#define _rotacion(w,u,v) \
     w.a0 = u.a0 * v.a0 - u.a1 * v.a1 \
            - u.a2 * v.a2 - u.a3 * v.a3;\
     w.a1 = u.a0 * v.a1 + u.a1 * v.a0 \
            + u.a2 * v.a3 - u.a3 * v.a2;\
     w.a2 = u.a0 * v.a2 - u.a1 * v.a3 \
            + u.a2 * v.a0 + u.a3 * v.a1;\
     w.a3 = u.a0 * v.a3 + u.a1 * v.a2 \
            - u.a2 * v.a1 + u.a3 * v.a0;


