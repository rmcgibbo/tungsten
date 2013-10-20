CC = mpic++
NETCDF_DIR = /Users/rmcgibbo/anaconda/pkgs/libnetcdf-4.2.1.1-1
OPENMM_LINKS = -L$(OPENMM_LIB_PATH) -lOpenMM -lOpenMMSerialization
NETCDF_LINKS = -L$(NETCDF_DIR)/lib -lnetcdf 

run: accelerator system.xml integrator.xml state.xml config.ini
	mpirun -np 2 ./accelerator system.xml integrator.xml state.xml config.ini

accelerator: build/mainloop.o build/inireader.o build/ini.o build/netcdfwriter.o
	$(CC) -o accelerator build/mainloop.o build/inireader.o build/ini.o $(OPENMM_LINKS) $(NETCDF_LINKS)

build/mainloop.o: src/mainloop.cpp
	$(CC) -o $@ -c $< -Iinclude -I$(OPENMM_INCLUDE_PATH) -Iinclude

build/ini.o: src/ini.c
	$(CC) -o $@ -c $< -Iinclude

build/inireader.o: src/INIReader.cpp
	$(CC) -o $@ -c $< -Iinclude

build/netcdfwriter.o: src/netcdfwriter.cpp
	$(CC) -o $@ -c $< -Iinclude -I$(NETCDF_DIR)/include

clean:
	rm -f build/* accelerator
