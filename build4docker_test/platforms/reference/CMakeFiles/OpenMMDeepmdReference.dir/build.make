# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/conda/bin/cmake

# The command to remove a file.
RM = /opt/conda/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/build4docker_test

# Include any dependencies generated for this target.
include platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/compiler_depend.make

# Include the progress variables for this target.
include platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/progress.make

# Include the compile flags for this target's objects.
include platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/flags.make

platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.o: platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/flags.make
platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.o: ../platforms/reference/src/ReferenceDeepmdKernelFactory.cpp
platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.o: platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/build4docker_test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.o"
	cd /mnt/build4docker_test/platforms/reference && /opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.o -MF CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.o.d -o CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.o -c /mnt/platforms/reference/src/ReferenceDeepmdKernelFactory.cpp

platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.i"
	cd /mnt/build4docker_test/platforms/reference && /opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/platforms/reference/src/ReferenceDeepmdKernelFactory.cpp > CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.i

platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.s"
	cd /mnt/build4docker_test/platforms/reference && /opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/platforms/reference/src/ReferenceDeepmdKernelFactory.cpp -o CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.s

platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.o: platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/flags.make
platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.o: ../platforms/reference/src/ReferenceDeepmdKernels.cpp
platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.o: platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/build4docker_test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.o"
	cd /mnt/build4docker_test/platforms/reference && /opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.o -MF CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.o.d -o CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.o -c /mnt/platforms/reference/src/ReferenceDeepmdKernels.cpp

platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.i"
	cd /mnt/build4docker_test/platforms/reference && /opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/platforms/reference/src/ReferenceDeepmdKernels.cpp > CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.i

platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.s"
	cd /mnt/build4docker_test/platforms/reference && /opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/platforms/reference/src/ReferenceDeepmdKernels.cpp -o CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.s

# Object files for target OpenMMDeepmdReference
OpenMMDeepmdReference_OBJECTS = \
"CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.o" \
"CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.o"

# External object files for target OpenMMDeepmdReference
OpenMMDeepmdReference_EXTERNAL_OBJECTS =

libOpenMMDeepmdReference.so: platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernelFactory.cpp.o
libOpenMMDeepmdReference.so: platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/src/ReferenceDeepmdKernels.cpp.o
libOpenMMDeepmdReference.so: platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/build.make
libOpenMMDeepmdReference.so: libOpenMMDeepmd.so
libOpenMMDeepmdReference.so: platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/build4docker_test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX shared library ../../libOpenMMDeepmdReference.so"
	cd /mnt/build4docker_test/platforms/reference && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/OpenMMDeepmdReference.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/build: libOpenMMDeepmdReference.so
.PHONY : platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/build

platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/clean:
	cd /mnt/build4docker_test/platforms/reference && $(CMAKE_COMMAND) -P CMakeFiles/OpenMMDeepmdReference.dir/cmake_clean.cmake
.PHONY : platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/clean

platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/depend:
	cd /mnt/build4docker_test && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt /mnt/platforms/reference /mnt/build4docker_test /mnt/build4docker_test/platforms/reference /mnt/build4docker_test/platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : platforms/reference/CMakeFiles/OpenMMDeepmdReference.dir/depend
