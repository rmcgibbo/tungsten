# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/rmcgibbo/projects/tungsten

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/rmcgibbo/projects/tungsten

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/rmcgibbo/projects/tungsten/CMakeFiles /home/rmcgibbo/projects/tungsten/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/rmcgibbo/projects/tungsten/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named tungsten

# Build rule for target.
tungsten: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tungsten
.PHONY : tungsten

# fast build rule for target.
tungsten/fast:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/build
.PHONY : tungsten/fast

src/INIReader.o: src/INIReader.cpp.o
.PHONY : src/INIReader.o

# target to build an object file
src/INIReader.cpp.o:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/INIReader.cpp.o
.PHONY : src/INIReader.cpp.o

src/INIReader.i: src/INIReader.cpp.i
.PHONY : src/INIReader.i

# target to preprocess a source file
src/INIReader.cpp.i:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/INIReader.cpp.i
.PHONY : src/INIReader.cpp.i

src/INIReader.s: src/INIReader.cpp.s
.PHONY : src/INIReader.s

# target to generate assembly for a file
src/INIReader.cpp.s:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/INIReader.cpp.s
.PHONY : src/INIReader.cpp.s

src/NetCDFTrajectoryFile.o: src/NetCDFTrajectoryFile.cpp.o
.PHONY : src/NetCDFTrajectoryFile.o

# target to build an object file
src/NetCDFTrajectoryFile.cpp.o:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/NetCDFTrajectoryFile.cpp.o
.PHONY : src/NetCDFTrajectoryFile.cpp.o

src/NetCDFTrajectoryFile.i: src/NetCDFTrajectoryFile.cpp.i
.PHONY : src/NetCDFTrajectoryFile.i

# target to preprocess a source file
src/NetCDFTrajectoryFile.cpp.i:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/NetCDFTrajectoryFile.cpp.i
.PHONY : src/NetCDFTrajectoryFile.cpp.i

src/NetCDFTrajectoryFile.s: src/NetCDFTrajectoryFile.cpp.s
.PHONY : src/NetCDFTrajectoryFile.s

# target to generate assembly for a file
src/NetCDFTrajectoryFile.cpp.s:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/NetCDFTrajectoryFile.cpp.s
.PHONY : src/NetCDFTrajectoryFile.cpp.s

src/ParallelKCenters.o: src/ParallelKCenters.cpp.o
.PHONY : src/ParallelKCenters.o

# target to build an object file
src/ParallelKCenters.cpp.o:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/ParallelKCenters.cpp.o
.PHONY : src/ParallelKCenters.cpp.o

src/ParallelKCenters.i: src/ParallelKCenters.cpp.i
.PHONY : src/ParallelKCenters.i

# target to preprocess a source file
src/ParallelKCenters.cpp.i:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/ParallelKCenters.cpp.i
.PHONY : src/ParallelKCenters.cpp.i

src/ParallelKCenters.s: src/ParallelKCenters.cpp.s
.PHONY : src/ParallelKCenters.s

# target to generate assembly for a file
src/ParallelKCenters.cpp.s:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/ParallelKCenters.cpp.s
.PHONY : src/ParallelKCenters.cpp.s

src/ParallelMSM.o: src/ParallelMSM.cpp.o
.PHONY : src/ParallelMSM.o

# target to build an object file
src/ParallelMSM.cpp.o:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/ParallelMSM.cpp.o
.PHONY : src/ParallelMSM.cpp.o

src/ParallelMSM.i: src/ParallelMSM.cpp.i
.PHONY : src/ParallelMSM.i

# target to preprocess a source file
src/ParallelMSM.cpp.i:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/ParallelMSM.cpp.i
.PHONY : src/ParallelMSM.cpp.i

src/ParallelMSM.s: src/ParallelMSM.cpp.s
.PHONY : src/ParallelMSM.s

# target to generate assembly for a file
src/ParallelMSM.cpp.s:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/ParallelMSM.cpp.s
.PHONY : src/ParallelMSM.cpp.s

src/csparse.o: src/csparse.c.o
.PHONY : src/csparse.o

# target to build an object file
src/csparse.c.o:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/csparse.c.o
.PHONY : src/csparse.c.o

src/csparse.i: src/csparse.c.i
.PHONY : src/csparse.i

# target to preprocess a source file
src/csparse.c.i:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/csparse.c.i
.PHONY : src/csparse.c.i

src/csparse.s: src/csparse.c.s
.PHONY : src/csparse.s

# target to generate assembly for a file
src/csparse.c.s:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/csparse.c.s
.PHONY : src/csparse.c.s

src/ini.o: src/ini.c.o
.PHONY : src/ini.o

# target to build an object file
src/ini.c.o:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/ini.c.o
.PHONY : src/ini.c.o

src/ini.i: src/ini.c.i
.PHONY : src/ini.i

# target to preprocess a source file
src/ini.c.i:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/ini.c.i
.PHONY : src/ini.c.i

src/ini.s: src/ini.c.s
.PHONY : src/ini.s

# target to generate assembly for a file
src/ini.c.s:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/ini.c.s
.PHONY : src/ini.c.s

src/mainloop.o: src/mainloop.cpp.o
.PHONY : src/mainloop.o

# target to build an object file
src/mainloop.cpp.o:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/mainloop.cpp.o
.PHONY : src/mainloop.cpp.o

src/mainloop.i: src/mainloop.cpp.i
.PHONY : src/mainloop.i

# target to preprocess a source file
src/mainloop.cpp.i:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/mainloop.cpp.i
.PHONY : src/mainloop.cpp.i

src/mainloop.s: src/mainloop.cpp.s
.PHONY : src/mainloop.s

# target to generate assembly for a file
src/mainloop.cpp.s:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/mainloop.cpp.s
.PHONY : src/mainloop.cpp.s

src/rmsd/theobald_rmsd.o: src/rmsd/theobald_rmsd.c.o
.PHONY : src/rmsd/theobald_rmsd.o

# target to build an object file
src/rmsd/theobald_rmsd.c.o:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/rmsd/theobald_rmsd.c.o
.PHONY : src/rmsd/theobald_rmsd.c.o

src/rmsd/theobald_rmsd.i: src/rmsd/theobald_rmsd.c.i
.PHONY : src/rmsd/theobald_rmsd.i

# target to preprocess a source file
src/rmsd/theobald_rmsd.c.i:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/rmsd/theobald_rmsd.c.i
.PHONY : src/rmsd/theobald_rmsd.c.i

src/rmsd/theobald_rmsd.s: src/rmsd/theobald_rmsd.c.s
.PHONY : src/rmsd/theobald_rmsd.s

# target to generate assembly for a file
src/rmsd/theobald_rmsd.c.s:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/rmsd/theobald_rmsd.c.s
.PHONY : src/rmsd/theobald_rmsd.c.s

src/utilities.o: src/utilities.cpp.o
.PHONY : src/utilities.o

# target to build an object file
src/utilities.cpp.o:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/utilities.cpp.o
.PHONY : src/utilities.cpp.o

src/utilities.i: src/utilities.cpp.i
.PHONY : src/utilities.i

# target to preprocess a source file
src/utilities.cpp.i:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/utilities.cpp.i
.PHONY : src/utilities.cpp.i

src/utilities.s: src/utilities.cpp.s
.PHONY : src/utilities.s

# target to generate assembly for a file
src/utilities.cpp.s:
	$(MAKE) -f CMakeFiles/tungsten.dir/build.make CMakeFiles/tungsten.dir/src/utilities.cpp.s
.PHONY : src/utilities.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... tungsten"
	@echo "... src/INIReader.o"
	@echo "... src/INIReader.i"
	@echo "... src/INIReader.s"
	@echo "... src/NetCDFTrajectoryFile.o"
	@echo "... src/NetCDFTrajectoryFile.i"
	@echo "... src/NetCDFTrajectoryFile.s"
	@echo "... src/ParallelKCenters.o"
	@echo "... src/ParallelKCenters.i"
	@echo "... src/ParallelKCenters.s"
	@echo "... src/ParallelMSM.o"
	@echo "... src/ParallelMSM.i"
	@echo "... src/ParallelMSM.s"
	@echo "... src/csparse.o"
	@echo "... src/csparse.i"
	@echo "... src/csparse.s"
	@echo "... src/ini.o"
	@echo "... src/ini.i"
	@echo "... src/ini.s"
	@echo "... src/mainloop.o"
	@echo "... src/mainloop.i"
	@echo "... src/mainloop.s"
	@echo "... src/rmsd/theobald_rmsd.o"
	@echo "... src/rmsd/theobald_rmsd.i"
	@echo "... src/rmsd/theobald_rmsd.s"
	@echo "... src/utilities.o"
	@echo "... src/utilities.i"
	@echo "... src/utilities.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

