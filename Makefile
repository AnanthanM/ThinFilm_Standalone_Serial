NAME 	:= EXECUTABLE
OBDIR	:= OBJ/
RUNDIR	:= RUN/

OBJECTS =  $(OBDIR)Thinfilm.o $(OBDIR)Cvode.o $(OBDIR)Sundials_nvector.o $(OBDIR)Sundials_math.o $(OBDIR)Cvode_spgmr.o $(OBDIR)Cvode_spils.o $(OBDIR)Sundials_spgmr.o $(OBDIR)Sundials_iterative.o $(OBDIR)Nvector_serial.o $(OBDIR)Sundials_band.o $(OBDIR)Sundials_direct.o $(OBDIR)Cvode_io.o

INCPATH	=	-I ~/ 

LIBPATH	= 	-L ~/

#LIBS	=	-lm /usr/lib/x86_64-linux-gnu/librt.so
LIBS	=	-lm 

HEADER 	:= CVODE.h SUNDIALS_NVECTOR.h SUNDIALS_TYPES.h SUNDIALS_CONFIG.h CVODE_IMPL.h SUNDIALS_MATH.h CVODE_SPGMR.h CVODE_SPILS.h SUNDIALS_SPGMR.h SUNDIALS_ITERATIVE.h CVODE_SPILS_IMPL.h NVECTOR_SERIAL.h SUNDIALS_BAND.h SUNDIALS_DIRECT.h      

CC	:=	gcc 

FLAGS	:=	$(INCPATH) -g  -Wall -O3 -march=native -ftree-vectorize -ftree-vectorizer-verbose=1 

LFLAGS	:=	$(LIBPATH) $(LIBS) 

####################################################################

$(NAME): $(OBJECTS) $(RUNDIR)
	$(CC) -o $(NAME) $(OBJECTS) $(FLAGS) $(LFLAGS)
	@mv $(NAME) $(RUNDIR) 

$(RUNDIR): 
	@test -d $(RUNDIR) || mkdir $(RUNDIR)

$(OBDIR)Thinfilm.o: Thinfilm.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Thinfilm.o Thinfilm.c

$(OBDIR)Cvode.o: Cvode.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Cvode.o Cvode.c

$(OBDIR)Sundials_nvector.o: Sundials_nvector.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Sundials_nvector.o Sundials_nvector.c

$(OBDIR)Sundials_math.o: Sundials_math.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Sundials_math.o Sundials_math.c

$(OBDIR)Cvode_spgmr.o: Cvode_spgmr.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Cvode_spgmr.o Cvode_spgmr.c

$(OBDIR)Cvode_spils.o: Cvode_spils.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Cvode_spils.o Cvode_spils.c

$(OBDIR)Sundials_spgmr.o: Sundials_spgmr.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Sundials_spgmr.o Sundials_spgmr.c

$(OBDIR)Sundials_iterative.o: Sundials_iterative.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Sundials_iterative.o Sundials_iterative.c

$(OBDIR)Nvector_serial.o: Nvector_serial.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Nvector_serial.o Nvector_serial.c

$(OBDIR)Sundials_band.o: Sundials_band.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Sundials_band.o Sundials_band.c

$(OBDIR)Sundials_direct.o: Sundials_direct.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Sundials_direct.o Sundials_direct.c

$(OBDIR)Cvode_io.o: Cvode_io.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Cvode_io.o Cvode_io.c

clean: 		
	rm -f *~
	rm $(OBDIR)*.o

Delete_Results: 		
	rm -f *~
	rm -r  $(RUNDIR)
	rm -r  $(OBDIR)
#	rm $(RUNDIR)/*.vtk
#	rm $(RUNDIR)/$(NAME)
