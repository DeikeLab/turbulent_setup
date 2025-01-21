#SET FILES OF INTEREST 
OUT_FOLDER = $(addprefix /home/nn8802/Documents/simulation/basilisk_dotc/mkf_,$(OUTPUT))
SRC_FILES = $(addsuffix .c,$(OUTPUT))
MPI_FILES = $(addprefix _,$(SRC_FILES))
OBJ_FILES = $(addsuffix .o,$(basename $(SRC_FILES)))

#SET C FLAGS AND BASILISK COMPILER, LINKING FLAGS
CC = gcc
CFLAGS += -events -source -Wall -O2 -D_MPI=1 -DDEBUG=1 
MPICCFLAGS += -Wall -std=c99 -O2 -D_MPI=1 -DDEBUG=1
BAS=qcc
MPICC = mpicc
#LDFLAGS = -lm -ldl
#LDFLAGS = -lgfortran -L${BASILISK}/ppr -lppr -lm -ldl
LDFLAGS = -lgfortran -L${BASILISK}/ppr -lm -ldl
OPENGLINC = ${BASILISK}
OPENGLIBS = -L${BASILISK}/gl -lglutils -lfb_osmesa -lGLU -lOSMesa


#MAKE COMMAND: TAKES ARGUMENT AND MAKES EXECUTABLE IN DESIRED FOLDER
$(OUTPUT): 
	$(BAS) $(CFLAGS) $(SRC_FILES)
	$(MPICC) $(MPICCFLAGS) $(MPI_FILES) -I$(OPENGLINC) -o $(OUT_FOLDER)/$(OUTPUT) $(OPENGLIBS) $(LDFLAGS)

