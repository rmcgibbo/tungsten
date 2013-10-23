CPP = mpic++
CC = mpicc

NETCDF_DIR = /usr
NETCDF_DIR = /Library/Frameworks/EPD64.framework/Versions/7.1
OPENMM_LINKS = -L$(OPENMM_LIB_PATH) -lOpenMM -lOpenMMSerialization
NETCDF_LINKS = -L$(NETCDF_DIR)/lib -lnetcdf_c++

INCLUDE = -Iinclude -Iinclude/rmsd -I$(OPENMM_INCLUDE_PATH) -I$(NETCDF_DIR)/include
CPP_FLAGS = $(INCLUDE) -fopenmp -msse2 -mssse3 -O0 -g
LD_FLAGS = $(OPENMM_LINKS) $(NETCDF_LINKS) -lgomp

CPP_FILES := $(wildcard src/*.cpp)
C_FILES := $(wildcard src/*.c) $(wildcard src/*/*.c)
CPP_OBJ_FILES = $(patsubst src/%.cpp,obj/%.o,$(CPP_FILES))
C_OBJ_FILES = $(patsubst src/%.c,obj/%.o,$(C_FILES))


tungsten: $(CPP_OBJ_FILES) $(C_OBJ_FILES)
	$(CPP) -o $@ $^ $(LD_FLAGS)

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CPP) $(CPP_FLAGS) -c -o $@ $<

obj/%.o: src/%.c
	@mkdir -p $(@D)
	$(CC) $(CPP_FLAGS) -c -o $@ $<

clean:
	rm -rf obj/* accelerator trj-*

run: tungsten
	mpirun -np 2 ./tungsten data/system.xml data/integrator.xml data/state.xml data/config.ini
