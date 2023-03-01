# Variables
CC = g++
VERSION = -std=c++23
OPT = -O3
DEBUG ?= no
MEMCHECK ?=

# Libraries settings
LIB_NAMES := Gillespie Monte_Carlo
LIB_SRC_DIRS := Shared_Objects
LIB_INCLUDE_DIRS := Shared_Objects $(LIB_NAMES)
LIB_SRCS := $(foreach dir, $(LIB_SRC_DIRS), $(wildcard $(dir)/*.cpp))
LIB_OBJS := $(foreach src, $(LIB_SRCS), $(patsubst %.cpp, %.o, $(src)))
LIB_HEADERS := $(foreach dir, $(LIB_INCLUDE_DIRS), $(wildcard $(dir)/*.h))
LIBS := $(foreach lib, $(LIB_NAMES), $(lib).so)

# Flags
ifeq ($(DEBUG), yes)
    FLAG = -DDEBUG
else
    FLAG =
endif

# Targets
all: $(LIBS)

define make_library
$(1).so: $(LIB_OBJS) $(1)/$(1).o
	$(CC) $(OPT) -shared -Wl,-soname,$(1).so -o $(1).so $(LIB_OBJS) $(1)/$(1).o $(FLAG)

$(1)/$(1).o: $(1)/$(1).cpp $(LIB_HEADERS)
	$(CC) $(VERSION) -fPIC $(OPT) -c $$< -o $$@ $(FLAG) $(MEMCHECK)
endef

$(foreach lib, $(LIB_NAMES), $(eval $(call make_library,$(lib))))

%.o: %.cpp $(LIB_HEADERS)
	$(CC) $(VERSION) -fPIC $(OPT) -c $< -o $@ $(FLAG) $(MEMCHECK)

.PHONY: clean EXEC

clean:
	rm -rf $(foreach dir, $(LIB_SRC_DIRS), $(dir)/*.o) *~ $(LIBS)

EXEC: $(LIB_OBJS)
	$(foreach lib, $(LIB_NAMES), $(CC) $(OPT) $(LIB_OBJS) $(lib)/$(lib).o -o $(lib).so $(FLAG) $(MEMCHECK);)
