# Flags
CFLAGS_1 = -Wall -I//Users/rodrigogonzalez/gsl/gsl-instal/include	# Compiling
CFLAGS_2 = -L/Users/rodrigogonzalez/gsl/gsl-instal/lib				# Linking
CFLAGS_3 = -lgsl -lgslcblas -lm										# Linking

all: executable 

# Linking
executable: main.o
	@echo "Linking and producing the final application"
	g++ $(CFLAGS_2) $(CFLAGS_3) main.o -o executable

# Compiling
main.o:
	@echo "Compiling the main file"
	g++ $(CFLAGS_1) -c ../src/main.cpp -o main.o

# Clear
clear:
	@echo "Removing all the files"
	@rm main.o executable
	 