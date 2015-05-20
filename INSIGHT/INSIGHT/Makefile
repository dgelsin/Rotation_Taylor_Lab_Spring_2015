############################################################
# Makefile for INSIGHT-EM
# - source code should be placed in  src/ subdirectory
# - also make sure to make subdirectories  obj/ and bin/
# - make all: generates latest version of INSIGHT-EM-v??? executable
# - make EM-leo: generates executable for previous version
# - make clean: cleans up all .o files and executables
############################################################

#####################
# COMPILATION SETUP #
#####################

# compiler
CC=gcc

# GSL - SET THE GSL DIRECTORY APPROPRIATELY
#
# standard setting (should work on swiftgen and local machines where GSL is installed)
GSLDIR=/usr/local/
#
# cac 
# GSLDIR=/cac/contrib/gsl-1.15/


# compilation options
#Debugging
#CFLAGS += -g -O0 -fstack-protector-all -Wall -DDEBUG
#Production
CFLAGS+= -fstack-protector-all -Wall -O3 -I$(GSLDIR)/include/
LINKFLAGS= -L$(GSLDIR)/lib/  $(GSLDIR)/lib/libgsl.a $(GSLDIR)/lib/libgsl.a -lm

ifeq ($(TARGETOS), Windows)
  CFLAGS += -DWINDOWS -liberty
endif

###########
# VERSION #
###########

version=1.1

########
# MAKE #
########

exec=INSIGHT-EM-v${version}

# for PRODUCTION
# --------------
all: bin/${exec}

# for DEBUG
# ---------
#all:
#	@echo 'Warning: Version $(version) is in development. Run "make debug" to complie, or upload previous version.'

#debug: obj/INSIGHT-EM.o obj/Utils.o obj/SumLogs.o obj/NumericOpt.o obj/bfgs.o
#	@echo -n 'Building debug version of executable bin/'${exec}'  -->  '
#	$(CC) obj/INSIGHT-EM.o obj/Utils.o obj/SumLogs.o obj/NumericOpt.o obj/bfgs.o $(CFLAGS) $(LINKFLAGS) -o bin/${exec}

clean:
	@echo -n 'Cleaning up exec and obj files            -->  '
	rm -rf obj/*.o bin/${exec}

bin/${exec}: obj/INSIGHT-EM.o obj/Utils.o obj/SumLogs.o obj/NumericOpt.o obj/bfgs.o
	@echo -n 'Building   executable bin/'${exec}' -->  '
	$(CC) obj/INSIGHT-EM.o obj/Utils.o obj/SumLogs.o obj/NumericOpt.o obj/bfgs.o $(CFLAGS) $(LINKFLAGS) -o bin/${exec}

obj/INSIGHT-EM.o: src/INSIGHT-EM.c src/Utils.h src/SumLogs.h 
	@echo -n 'Compliling source file src/INSIGHT-EM.c   -->  '
	$(CC) $(CFLAGS) -c src/INSIGHT-EM.c -o obj/INSIGHT-EM.o

obj/Utils.o: src/Utils.c src/Utils.h
	@echo -n 'Compliling source file src/Utils.c        -->  '
	$(CC) $(CFLAGS) -c src/Utils.c -o obj/Utils.o

obj/SumLogs.o: src/SumLogs.c src/SumLogs.h src/NumericOpt.h 
	@echo -n 'Compliling source file src/SumLogs.c      -->  '
	$(CC) $(CFLAGS) -c src/SumLogs.c -o obj/SumLogs.o

obj/NumericOpt.o: src/NumericOpt.c src/NumericOpt.h src/bfgs.h 
	@echo -n 'Compliling source file src/NumericOpt.c   -->  '
	$(CC) $(CFLAGS) -c src/NumericOpt.c -o obj/NumericOpt.o

obj/bfgs.o: src/bfgs.c src/bfgs.h 
	@echo -n 'Compliling source file src/bfgs.c         -->  '
	$(CC) $(CFLAGS) -c src/bfgs.c -o obj/bfgs.o

