# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.3

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /media/hoquangnam/Data/Xilinx/SDK/2018.3/tps/lnx64/cmake-3.3.2/bin/cmake

# The command to remove a file.
RM = /media/hoquangnam/Data/Xilinx/SDK/2018.3/tps/lnx64/cmake-3.3.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/hoquangnam/Documents/Triangulation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/hoquangnam/Documents/Triangulation/build

# Include any dependencies generated for this target.
include CMakeFiles/triangulation.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/triangulation.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/triangulation.dir/flags.make

CMakeFiles/triangulation.dir/triangulation.o: CMakeFiles/triangulation.dir/flags.make
CMakeFiles/triangulation.dir/triangulation.o: ../triangulation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/hoquangnam/Documents/Triangulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/triangulation.dir/triangulation.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/triangulation.dir/triangulation.o -c /home/hoquangnam/Documents/Triangulation/triangulation.cpp

CMakeFiles/triangulation.dir/triangulation.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/triangulation.dir/triangulation.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/hoquangnam/Documents/Triangulation/triangulation.cpp > CMakeFiles/triangulation.dir/triangulation.i

CMakeFiles/triangulation.dir/triangulation.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/triangulation.dir/triangulation.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/hoquangnam/Documents/Triangulation/triangulation.cpp -o CMakeFiles/triangulation.dir/triangulation.s

CMakeFiles/triangulation.dir/triangulation.o.requires:

.PHONY : CMakeFiles/triangulation.dir/triangulation.o.requires

CMakeFiles/triangulation.dir/triangulation.o.provides: CMakeFiles/triangulation.dir/triangulation.o.requires
	$(MAKE) -f CMakeFiles/triangulation.dir/build.make CMakeFiles/triangulation.dir/triangulation.o.provides.build
.PHONY : CMakeFiles/triangulation.dir/triangulation.o.provides

CMakeFiles/triangulation.dir/triangulation.o.provides.build: CMakeFiles/triangulation.dir/triangulation.o


CMakeFiles/triangulation.dir/Libs/coordinate.o: CMakeFiles/triangulation.dir/flags.make
CMakeFiles/triangulation.dir/Libs/coordinate.o: ../Libs/coordinate.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/hoquangnam/Documents/Triangulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/triangulation.dir/Libs/coordinate.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/triangulation.dir/Libs/coordinate.o -c /home/hoquangnam/Documents/Triangulation/Libs/coordinate.cpp

CMakeFiles/triangulation.dir/Libs/coordinate.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/triangulation.dir/Libs/coordinate.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/hoquangnam/Documents/Triangulation/Libs/coordinate.cpp > CMakeFiles/triangulation.dir/Libs/coordinate.i

CMakeFiles/triangulation.dir/Libs/coordinate.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/triangulation.dir/Libs/coordinate.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/hoquangnam/Documents/Triangulation/Libs/coordinate.cpp -o CMakeFiles/triangulation.dir/Libs/coordinate.s

CMakeFiles/triangulation.dir/Libs/coordinate.o.requires:

.PHONY : CMakeFiles/triangulation.dir/Libs/coordinate.o.requires

CMakeFiles/triangulation.dir/Libs/coordinate.o.provides: CMakeFiles/triangulation.dir/Libs/coordinate.o.requires
	$(MAKE) -f CMakeFiles/triangulation.dir/build.make CMakeFiles/triangulation.dir/Libs/coordinate.o.provides.build
.PHONY : CMakeFiles/triangulation.dir/Libs/coordinate.o.provides

CMakeFiles/triangulation.dir/Libs/coordinate.o.provides.build: CMakeFiles/triangulation.dir/Libs/coordinate.o


# Object files for target triangulation
triangulation_OBJECTS = \
"CMakeFiles/triangulation.dir/triangulation.o" \
"CMakeFiles/triangulation.dir/Libs/coordinate.o"

# External object files for target triangulation
triangulation_EXTERNAL_OBJECTS =

triangulation: CMakeFiles/triangulation.dir/triangulation.o
triangulation: CMakeFiles/triangulation.dir/Libs/coordinate.o
triangulation: CMakeFiles/triangulation.dir/build.make
triangulation: CMakeFiles/triangulation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/hoquangnam/Documents/Triangulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable triangulation"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/triangulation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/triangulation.dir/build: triangulation

.PHONY : CMakeFiles/triangulation.dir/build

CMakeFiles/triangulation.dir/requires: CMakeFiles/triangulation.dir/triangulation.o.requires
CMakeFiles/triangulation.dir/requires: CMakeFiles/triangulation.dir/Libs/coordinate.o.requires

.PHONY : CMakeFiles/triangulation.dir/requires

CMakeFiles/triangulation.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/triangulation.dir/cmake_clean.cmake
.PHONY : CMakeFiles/triangulation.dir/clean

CMakeFiles/triangulation.dir/depend:
	cd /home/hoquangnam/Documents/Triangulation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/hoquangnam/Documents/Triangulation /home/hoquangnam/Documents/Triangulation /home/hoquangnam/Documents/Triangulation/build /home/hoquangnam/Documents/Triangulation/build /home/hoquangnam/Documents/Triangulation/build/CMakeFiles/triangulation.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/triangulation.dir/depend

