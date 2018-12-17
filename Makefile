#CFLAGS= -g -Wall                             #para debug  
CFLAGS=  -O2 -funroll-loops -finline-functions #optimizacion maxima

OBJ   = su2h_ini.o su2h_main.o su2h_med.o su2h_upd.o medidas.o

su2h: $(OBJ)
	gcc $(CFLAGS) $(OBJ) -lm  -o Su2H
.c.o:
	gcc -DSUN4 -c $(CFLAGS) $<

clean: 
	/bin/rm -f $(OBJ) Su2H

#		*Individual File Dependencies*
su2h_main.o: su2h_main.c su2h.h su2h_arit.h 

su2h_ini.o: su2h_ini.c  su2h.h su2h_arit.h

su2h_med.o: su2h_med.c su2h.h su2h_arit.h

su2h_upd.o: su2h_upd.c su2h.h su2h_arit.h

medidas.o: medidas.c su2h.h su2h_arit.h

