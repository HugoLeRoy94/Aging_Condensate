CC =g++ #compilateur
VERSION = -std=c++23#derniere version de c++
OPT = -O3 #optimisation de ouf
SharedLib = lib.so
DEBUG=yes
NAME=Aight
MEMCHECK= #-g

SRC=$(wildcard *.cpp)# $(wildcard *.c)


OBJ=$(SRC:.cpp=.o) #convertie la liste de .cpp en une liste de .o

HEAD=$(filter-out main.h System_c.h $(wildcard *.c) , $(SRC:.cpp=.h)) #converti tout les .cpp en .h sauf main.cpp et certains autres...

ifeq ($(DEBUG),yes)
	FLAG=-DDEBUG
else
	FLAG=
endif

all : $(SharedLib)
$(SharedLib) : $(OBJ)
	     $(CC) $(OPT) -shared -Wl,-soname,$(SharedLib) -o $(SharedLib) $(OBJ) $(FLAG)
main.o: $(SRC) $(HEAD)

.cpp.o:
	$(CC) $(VERSION) -fPIC $(OPT) -c $< -o $@ $(FLAG) $(MEMCHECK)

.PHONY: clean

clean:
	rm -rf *.o *~ $(SharedLib)
.PHONY: EXEC

EXEC: $(OBJ)
	$(CC) $(OPT) $(OBJ)  -o $(NAME) $(FLAG) $(MEMCHECK)
