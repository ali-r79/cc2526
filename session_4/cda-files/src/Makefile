#FC=g95
#STDF= -std=F
FC=gfortran
FFLAGS =

TARGET= ../bin/cda 
OBJECTS= cda.o \
         kinds.o \
         cubes.o

.SUFFIXES:
.SUFFIXES: .f95 .o

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJECTS)

cda.o: cda.f95 kinds.o cubes.o
cubes.o: cubes.f95 kinds.o
kinds.o: kinds.f95

%.o: %.f95
	$(FC) $(STDF) $(FFLAGS) -c -I. -o $@ $<

.PHONY: clean
clean:
	$(RM) $(TARGET) *.mod *.o
