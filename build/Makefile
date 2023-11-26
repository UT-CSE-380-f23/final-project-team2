# Flags
#CFLAGS_1 = -Wall -I /home/kenneth/soft/gsl-install/include # Compiling
#CFLAGS_2 = -L /home/kenneth/soft/gsl-install/lib # Linking
#CFLAGS_1 = -Wall -I$(WORK)/soft/bin/grvy/include -I$(TACC_GSL_INC) # Compiling
#CFLAGS_2 = -L $(WORK)/soft/bin/grvy/lib -L$(TACC_GSL_LIB) # Linking
#CFLAGS_3 = -lm -lgsl -lgslcblas -lgrvy # Linking

CFLAGS_1 = -Wall -I$(WORK)/soft/bin/grvy/include -I$(TACC_GSL_INC) # Compiling
CFLAGS_2 = -L $(WORK)/soft/bin/grvy/lib -L$(TACC_GSL_LIB) # Linking

#all: executable 

# Linking
#executable: main.o
#	@echo "Linking and producing the final application"
#	g++ $(CFLAGS_2) $(CFLAGS_3) main.o -o executable

# Compiling
#main.o:
#	@echo "Compiling the main file"
#	g++ $(CFLAGS_1) -c ../src/main.cpp -o main.o

# Clear
#clear:
#	@echo "Removing all the files"
#	@rm main.o executable
# Files
EXEC := main
SRC  := $(wildcard *.cpp)
OBJ  := $(patsubst %.cpp,%.o,$(SRC))

# Options
CC      := g++
CFLAGS  := -std=c++17 -O3
#LDFLAGS := -L /home/kenneth/soft/gsl-install/lib -I /home/kenneth/soft/gsl-install/include
#LDLIBS  := -lm -lgsl -lgslcblas
LDFLAGS := $(CFLAGS_1) $(CFLAGS_2) 
LDLIBS := -lgsl -lgslcblas -lm -lgrvy

# Rules
$(EXEC): $(OBJ)
	$(CC) $(LDFLAGS) -o $@  $^ $(LDLIBS)
%.o:  %.cpp
	$(CC) $(CFLAGS) $(LDFLAGS)  -c $< $(LDLIBS)
oned_solver.o : fd_solver.h

# Useful phony targets
.PHONY: clobber clean neat echo
clobber:  clean
	$(RM) $(EXEC)
clean:  neat
	$(RM) $(OBJ)
neat:
	$(RM) *~ .*~
echo:
	@echo $(OBJ)