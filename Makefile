# Get Trilinos as one entity

include /opt/local/Trilinos/include/Makefile.export.Trilinos


INCP	= -I/opt/local/include/newmat
LSPECTRAL	= -L/opt/local/lib -lnewmat
MY_RPATH = -rpath $(subst -L, ,$(Trilinos_LIBRARY_DIRS))

# Make sure to use same compilers and flags as Trilinos
CXX=$(Trilinos_CXX_COMPILER)
CC=$(Trilinos_C_COMPILER)
FORT=$(Trilinos_Fortran_COMPILER)

CXX_FLAGS=$(Trilinos_CXX_COMPILER_FLAGS) $(USER_CXX_FLAGS)
#-Wc++11-extensions
C_FLAGS=$(Trilinos_C_COMPILER_FLAGS) $(USERC_FLAGS)
FORT_FLAGS=$(Trilinos_Fortran_COMPILER_FLAGS) $(USER_FORT_FLAGS)

INCLUDE_DIRS=$(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS) $(INCP)
LIBRARY_DIRS=-L/usr/lib $(Trilinos_LIBRARY_DIRS)  $(Trilinos_TPL_LIBRARY_DIRS)
LIBRARIES=$(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) $(LSPECTRAL)

LINK_FLAGS=$(Trilinos_EXTRA_LD_FLAGS)

#just assuming that epetra is turned on.
#DEFINES=-DMYAPP_EPETRA

CLINKER	= $(CXX) 
##$(Trilinos_LINKER)


####### Target

TARGET	= main
DESTDIR = ./
VER_MAJ = 1
VER_MIN = 0

####### Arquivos que precisam ser mudados de acordo com o seu objetivo

TEMPLATED = 	GeProb.hpp \
                PhElem.hpp 	

HEADERS =	DG_EI_Header.h\
		DG_Prob.h\
		DG_Elem.hpp \
		DG_saidas_intermediarias.h\
		Fluids.h\
		Funcoes_c.h\
		Geo.h\
		Hexahedral.h\
		Linear.h\
		MyOptions.h\
		MyTrilinos.h\
		Particao.h\
		Quadrilateral.h\
		Stdel.h\
		Tetrahedral.h\
		Triangle.h\
		Tstruct.h\
		spectral.h
SOURCES =	AMPFunctionsA.cpp\
		ASPFunctions.cpp\
		DG.cpp\
		DG_EI_Inflow.cpp\
		DG_EI_Inflow_b.cpp\
		DG_EI_Inflow_a.cpp\
		DG_EI_Interior.cpp\
		DG_EI_Interior_b.cpp\
		DG_EI_Interior_a.cpp\
		DG_EI_Outflow.cpp\
		DG_EI_Outflow_b.cpp\
		DG_EI_Outflow_a.cpp\
		DG_EI_flux.cpp\
		DG_IG_EI.cpp\
		DG_Iterate.cpp\
		DG_MVRA.cpp\
		DG_Prob.cpp\
		DG_Elem.cpp \
		DG_driver.cpp\
		DG_eco.cpp\
		DG_preamble.cpp\
		DG_saidas_intermediarias.cpp\
		Differentiation.cpp\
		DG_Eigenvectors.cpp\
		Fluids.cpp\
		Funcoes_globais.cpp\
		Hexahedral.cpp\
		Linear.cpp\
		Particao.cpp\
		Quadrilateral.cpp\
		Stdel.cpp\
		Tetrahedral.cpp\
		Transfere_rst.cpp\
		Triangle.cpp\
		main.cpp
OBJECTS =	AMPFunctionsA.o\
		ASPFunctions.o\
		DG.o\
		DG_EI_Inflow.o\
		DG_EI_Inflow_b.o\
		DG_EI_Inflow_a.o\
		DG_EI_Interior.o\
		DG_EI_Interior_b.o\
		DG_EI_Interior_a.o\
		DG_EI_Outflow.o\
		DG_EI_Outflow_b.o\
		DG_EI_Outflow_a.o\
		DG_EI_flux.o\
		DG_IG_EI.o\
		DG_Iterate.o\
		DG_MVRA.o\
		DG_Prob.o\
		DG_Elem.o \
		DG_driver.o\
		DG_eco.o\
		DG_preamble.o\
		DG_saidas_intermediarias.o\
		Differentiation.o\
		DG_Eigenvectors.o\
		Fluids.o\
		Funcoes_globais.o\
		Hexahedral.o\
		Linear.o\
		Particao.o\
		Quadrilateral.o\
		Stdel.o\
		Tetrahedral.o\
		Transfere_rst.o\
		Triangle.o
####### Implicit rules

.SUFFIXES: .cpp .cxx .cc .C .c

.cpp.o:
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE_DIRS) -o $@ $<

.cxx.o:
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE_DIRS) -o $@ $<

.cc.o:
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE_DIRS) -o $@ $<

.C.o:
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE_DIRS) -o $@ $<

.c.o:
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE_DIRS) -o $@ $<

####### Build rules

all: $(DESTDIR)$(TARGET)


####### Linking

$(DESTDIR)$(TARGET): main.o $(OBJECTS) 
	$(CLINKER)    -o $(TARGET) main.o $(MY_RPATH) $(OBJECTS) $(LINK_FLAGS) $(LIBRARY_DIRS) $(LIBRARIES)


## A linha abaixo usa install_name_tool para adicionar o rpath ao executavel main	
#install_name_tool -add_rpath $(MY_RPATH) main

interpolar: ValoresInterpolados.o $(OBJECTS)
	    $(CLINKER)    -o interpolar ValoresInterpolados.o $(MY_RPATH) $(OBJECTS) $(LINK_FLAGS) $(LIBRARY_DIRS) $(LIBRARIES)

teste_ordenar: teste_ordenar.o $(OBJECTS)
	    $(CLINKER)    -o teste_ordenar teste_ordenar.o $(MY_RPATH) $(OBJECTS) $(LINK_FLAGS) $(LIBRARY_DIRS) $(LIBRARIES)


.PHONY: clean
clean:
	rm -f *.o *.a *.exe $(TARGET) *~ core

####### Compile 

$(OBJECTS): $(HEADERS) $(TEMPLATED)

$(OBJECTS_IP): $(HEADERS) $(TEMPLATED)

$(OBJECTS_CF): $(HEADERS) $(TEMPLATED)

main.o: main.cpp $(HEADERS) $(TEMPLATED)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE_DIRS) -o $@  main.cpp

ValoresInterpolados.o: ValoresInterpolados.cc $(HEADERS) $(TEMPLATED)
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE_DIRS) -o $@  ValoresInterpolados.cc

teste_ordenar.o: teste_ordenar.cc
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE_DIRS) -o $@  teste_ordenar.cc
