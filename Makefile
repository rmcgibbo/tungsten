CC = mpic++
NETCDF_DIR = /usr
OPENMM_LINKS = -L$(OPENMM_LIB_PATH) -lOpenMM -lOpenMMSerialization
NETCDF_LINKS = -L$(NETCDF_DIR)/lib -lnetcdf_c++

run: accelerator system.xml integrator.xml state.xml config.ini
	mpirun -np 2 ./accelerator system.xml integrator.xml state.xml config.ini

accelerator: build/mainloop.o build/INIReader.o build/ini.o build/NetCDFTrajectoryFile.o build/ParallelKCenters.o build/utilities.o
	$(CC) -o accelerator build/mainloop.o build/INIReader.o build/ini.o build/NetCDFTrajectoryFile.o build/ParallelKCenters.o build/utilities.o $(OPENMM_LINKS) $(NETCDF_LINKS)

build/mainloop.o: src/mainloop.cpp
	$(CC) -o $@ -c $< -Iinclude -I$(OPENMM_INCLUDE_PATH)

build/ini.o: src/ini.c
	$(CC) -o $@ -c $< -Iinclude

build/INIReader.o: src/INIReader.cpp
	$(CC) -o $@ -c $< -Iinclude

build/NetCDFTrajectoryFile.o: src/NetCDFTrajectoryFile.cpp
	$(CC) -o $@ -c $< -Iinclude -I$(NETCDF_DIR)/include -I$(OPENMM_INCLUDE_PATH)

build/ParallelKCenters.o: src/ParallelKCenters.cpp
	$(CC) -o $@ -c $< -Iinclude -I$(OPENMM_INCLUDE_PATH)


build/utilities.o: src/utilities.cpp
	$(CC) -o $@ -c $< -Iinclude -I$(OPENMM_INCLUDE_PATH)


clean:
	rm -rf build/* accelerator trj-*
