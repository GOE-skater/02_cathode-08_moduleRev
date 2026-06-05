# ----------------------------------------------------------------
# environment
CC = mpicxx

# ----------------------------------------------------------------
# options
CFLAGS = -std=c++17
CINCLUDES = -I ./include -I /opt/homebrew/Cellar/hypre/2.31.0/include -I /opt/homebrew/Cellar/open-mpi/5.0.3_1/include -I /opt/homebrew/include
CLIBS = -L /opt/homebrew/Cellar/hypre/2.31.0/lib -L/opt/homebrew/lib 
CLIBOPTIONS = -lHYPRE -lyaml-cpp
OPTIMIZE = -Ofast#-fopenmp
EXE_NAME = a.out

# ----------------------------------------------------------------
# executables
Normal: ./src/main.cpp
	$(CC) $(CFLAGS) $(CINCLUDES) $(CLIBS) $(OPTIMIZE) ./src/main.cpp -o $(EXE_NAME) $(CLIBOPTIONS)
