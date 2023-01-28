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
include CMakeFiles/OpenMMDeepmd.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/OpenMMDeepmd.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/OpenMMDeepmd.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/OpenMMDeepmd.dir/flags.make

CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.o: CMakeFiles/OpenMMDeepmd.dir/flags.make
CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.o: ../openmmapi/src/DeepmdForce.cpp
CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.o: CMakeFiles/OpenMMDeepmd.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/build4docker_test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.o"
	/opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.o -MF CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.o.d -o CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.o -c /mnt/openmmapi/src/DeepmdForce.cpp

CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.i"
	/opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/openmmapi/src/DeepmdForce.cpp > CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.i

CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.s"
	/opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/openmmapi/src/DeepmdForce.cpp -o CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.s

CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.o: CMakeFiles/OpenMMDeepmd.dir/flags.make
CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.o: ../openmmapi/src/DeepmdForceImpl.cpp
CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.o: CMakeFiles/OpenMMDeepmd.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/build4docker_test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.o"
	/opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.o -MF CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.o.d -o CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.o -c /mnt/openmmapi/src/DeepmdForceImpl.cpp

CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.i"
	/opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/openmmapi/src/DeepmdForceImpl.cpp > CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.i

CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.s"
	/opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/openmmapi/src/DeepmdForceImpl.cpp -o CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.s

CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.o: CMakeFiles/OpenMMDeepmd.dir/flags.make
CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.o: ../serialization/src/DeepmdForceProxy.cpp
CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.o: CMakeFiles/OpenMMDeepmd.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/build4docker_test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.o"
	/opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.o -MF CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.o.d -o CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.o -c /mnt/serialization/src/DeepmdForceProxy.cpp

CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.i"
	/opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/serialization/src/DeepmdForceProxy.cpp > CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.i

CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.s"
	/opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/serialization/src/DeepmdForceProxy.cpp -o CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.s

CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.o: CMakeFiles/OpenMMDeepmd.dir/flags.make
CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.o: ../serialization/src/DeepmdSerializationProxyRegistration.cpp
CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.o: CMakeFiles/OpenMMDeepmd.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/build4docker_test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.o"
	/opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.o -MF CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.o.d -o CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.o -c /mnt/serialization/src/DeepmdSerializationProxyRegistration.cpp

CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.i"
	/opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/serialization/src/DeepmdSerializationProxyRegistration.cpp > CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.i

CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.s"
	/opt/conda/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/serialization/src/DeepmdSerializationProxyRegistration.cpp -o CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.s

# Object files for target OpenMMDeepmd
OpenMMDeepmd_OBJECTS = \
"CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.o" \
"CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.o" \
"CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.o" \
"CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.o"

# External object files for target OpenMMDeepmd
OpenMMDeepmd_EXTERNAL_OBJECTS =

libOpenMMDeepmd.so: CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForce.cpp.o
libOpenMMDeepmd.so: CMakeFiles/OpenMMDeepmd.dir/openmmapi/src/DeepmdForceImpl.cpp.o
libOpenMMDeepmd.so: CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdForceProxy.cpp.o
libOpenMMDeepmd.so: CMakeFiles/OpenMMDeepmd.dir/serialization/src/DeepmdSerializationProxyRegistration.cpp.o
libOpenMMDeepmd.so: CMakeFiles/OpenMMDeepmd.dir/build.make
libOpenMMDeepmd.so: CMakeFiles/OpenMMDeepmd.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/build4docker_test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX shared library libOpenMMDeepmd.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/OpenMMDeepmd.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/OpenMMDeepmd.dir/build: libOpenMMDeepmd.so
.PHONY : CMakeFiles/OpenMMDeepmd.dir/build

CMakeFiles/OpenMMDeepmd.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/OpenMMDeepmd.dir/cmake_clean.cmake
.PHONY : CMakeFiles/OpenMMDeepmd.dir/clean

CMakeFiles/OpenMMDeepmd.dir/depend:
	cd /mnt/build4docker_test && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt /mnt /mnt/build4docker_test /mnt/build4docker_test /mnt/build4docker_test/CMakeFiles/OpenMMDeepmd.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/OpenMMDeepmd.dir/depend

