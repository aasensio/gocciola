# Makefile for building Gocciola

# Directories with libraries
BACKOPA_LIB = SOPAS

# Select compiler
# Intel Fortran compiler
COMPILER = ifort
OPTS = -c -O3 -warn all
#-check all -g -traceback -stand f90
OPTS_LINK = -O3
#LIBS = -L$(BACKOPA_LIB) -lbackopa

# GFortran
#COMPILER = gfortran
#OPTS = -c -O3 -Wall
#OPTS_LINK = -O3
# LIBS = -L$(BACKOPA_LIB) -lbackopa

OPTIONS = $(OPTS) $(DEBUG)

OBJECTS = accel.o background_opa.o dealloc.o formal_jac_overl.o math.o pp_jac_overl.o\
atmosfer.o emerging.o general.o overlap.o resultados.o\
atomic.o init.o pp_cep.o short_ch.o\
background.o formal_gs.o pp_gs.o statis.o\
collision.o formal_gs_overl.o lvg.o pp_gs_overl.o vars.o\
convol.o formal_jac.o pp_jac.o zeeman.o jac_molec.o 

	
gocciola: $(OBJECTS)
	$(COMPILER) $(OBJECTS) $(OPTS_LINK) -o gocciola $(LIBS)
	cp gocciola ../

all:
#	make sopas
	make gocciola 
	
# sopas:
# 	cd SOPAS; $(COMPILER) $(OPTS) *.f ; ar cr libbackopa.a *.o

clean:
	find . -maxdepth 4 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete
	find . -maxdepth 1 -name "*.f90~" -delete ; find . -maxdepth 1 -name "gocciola" -delete
	find ../ -maxdepth 1 -name "gocciola" -delete

accel.o : accel.f90 math.o statis.o vars.o
	$(COMPILER) $(OPTIONS) accel.f90
	
atmosfer.o : atmosfer.f90 general.o vars.o 
	$(COMPILER) $(OPTIONS) atmosfer.f90
	
background_opa.o : background_opa.f90
	$(COMPILER) $(OPTIONS) background_opa.f90
	
atomic.o : atomic.f90 collision.o general.o vars.o 
	$(COMPILER) $(OPTIONS) atomic.f90
	
background.o : background.f90 vars.o background_opa.o
	$(COMPILER) $(OPTIONS) background.f90
	
collision.o : collision.f90 vars.o 
	$(COMPILER) $(OPTIONS) collision.f90
	
convol.o : convol.f90 vars.o emerging.o
	$(COMPILER) $(OPTIONS) convol.f90
	
dealloc.o : dealloc.f90 vars.o 
	$(COMPILER) $(OPTIONS) dealloc.f90
	
emerging.o : emerging.f90 general.o pp_jac_overl.o pp_jac.o formal_jac_overl.o formal_jac.o short_ch.o vars.o 
	$(COMPILER) $(OPTIONS) emerging.f90
	
# flux.o : flux.f90 vars.o 
# 	$(COMPILER) $(OPTIONS) flux.f90
	
formal_gs.o : formal_gs.f90 short_ch.o statis.o general.o vars.o 
	$(COMPILER) $(OPTIONS) formal_gs.f90
	
formal_gs_overl.o : formal_gs_overl.f90 overlap.o short_ch.o statis.o general.o vars.o 
	$(COMPILER) $(OPTIONS) formal_gs_overl.f90
	
formal_jac.o : formal_jac.f90 short_ch.o general.o vars.o 
	$(COMPILER) $(OPTIONS) formal_jac.f90
	
formal_jac_overl.o : formal_jac_overl.f90 overlap.o short_ch.o general.o vars.o 
	$(COMPILER) $(OPTIONS) formal_jac_overl.f90
	
general.o : general.f90 math.o vars.o 
	$(COMPILER) $(OPTIONS) general.f90
	
init.o : init.f90 math.o resultados.o overlap.o background.o general.o atmosfer.o atomic.o vars.o 
	$(COMPILER) $(OPTIONS) init.f90
	
jac_molec.o : jac_molec.f90 zeeman.o emerging.o dealloc.o resultados.o vars.o vars.o accel.o atmosfer.o statis.o pp_gs_overl.o pp_gs.o formal_gs_overl.o formal_gs.o pp_cep.o pp_jac_overl.o pp_jac.o formal_jac_overl.o formal_jac.o init.o general.o 
	$(COMPILER) $(OPTIONS) jac_molec.f90
	
lvg.o : lvg.f90 math.o resultados.o atmosfer.o statis.o vars.o
	$(COMPILER) $(OPTIONS) lvg.f90
	
math.o : math.f90 vars.o 
	$(COMPILER) $(OPTIONS) math.f90
	
overlap.o : overlap.f90 general.o vars.o 
	$(COMPILER) $(OPTIONS) overlap.f90
	
pp_cep.o : pp_cep.f90 statis.o atmosfer.o math.o short_ch.o general.o vars.o 
	$(COMPILER) $(OPTIONS) pp_cep.f90
	
pp_gs.o : pp_gs.f90 short_ch.o statis.o general.o vars.o 
	$(COMPILER) $(OPTIONS) pp_gs.f90
	
pp_gs_overl.o : pp_gs_overl.f90 overlap.o short_ch.o statis.o general.o vars.o 
	$(COMPILER) $(OPTIONS) pp_gs_overl.f90
	
pp_jac.o : pp_jac.f90 math.o short_ch.o general.o vars.o 
	$(COMPILER) $(OPTIONS) pp_jac.f90
	
pp_jac_overl.o : pp_jac_overl.f90 overlap.o short_ch.o general.o vars.o 
	$(COMPILER) $(OPTIONS) pp_jac_overl.f90
	
resultados.o : resultados.f90 general.o emerging.o vars.o 
	$(COMPILER) $(OPTIONS) resultados.f90
	
short_ch.o : short_ch.f90 vars.o 
	$(COMPILER) $(OPTIONS) short_ch.f90
	
statis.o : statis.f90 math.o vars.o 
	$(COMPILER) $(OPTIONS) statis.f90
	
vars.o : vars.f90 
	$(COMPILER) $(OPTIONS) vars.f90
	
zeeman.o : zeeman.f90 general.o vars.o
	$(COMPILER) $(OPTIONS) zeeman.f90
	
