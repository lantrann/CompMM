.SUFFIXES: .f .F .o .a  .f90 .F90 .f95 .F95

default: what excbias

what:
	@ echo 
	@ echo 'Compiling the program for core-shell magnetic nanoparticles'
	@ echo 'The program is writen by L. B. Ho and T. N. Lan '
	@ echo 'Hit ^C to exit'
	@ echo 
	@ sleep 2

FC = gfortran
FCFLAGS =
LDFAGS =

SRCDIR=./Src  
VPATH=$(SRCDIR)/ExchangeBias $(SRCDIR)/Hyperthemia $(SRCDIR)/Nanodisks

OBJDIR=

EBOBJS = eb_module.o \
	 input.o \
	 set_info_site.o \
         small_trial_step.o \
         uniform_trial_step.o \
         reflect_trial_step.o \
         eta_max.o \
         metropolis.o \
	 uniaxial_anisotropy_z.o \
	 zeeman_z.o \
	 exchange_coupling.o \
	 set_init_config.o \
	 sweep_core.o \
	 sweep_shell.o \
	 hysteresis.o \
	 three_trial_steps.o \
	 eb_main.o

excbias: $(EBOBJS)
	 $(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f95
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean

clean:
	rm -f excbias
	rm -f *.o *.mod *.MOD
