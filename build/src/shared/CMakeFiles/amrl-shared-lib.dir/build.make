# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/skm2019/Documents/CS393R/cs393r_starter

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/skm2019/Documents/CS393R/cs393r_starter/build

# Include any dependencies generated for this target.
include src/shared/CMakeFiles/amrl-shared-lib.dir/depend.make

# Include the progress variables for this target.
include src/shared/CMakeFiles/amrl-shared-lib.dir/progress.make

# Include the compile flags for this target's objects.
include src/shared/CMakeFiles/amrl-shared-lib.dir/flags.make

src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o: src/shared/CMakeFiles/amrl-shared-lib.dir/flags.make
src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o: ../src/shared/util/helpers.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skm2019/Documents/CS393R/cs393r_starter/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o -c /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/helpers.cc

src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.i"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/helpers.cc > CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.i

src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.s"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/helpers.cc -o CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.s

src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o.requires:

.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o.requires

src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o.provides: src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o.requires
	$(MAKE) -f src/shared/CMakeFiles/amrl-shared-lib.dir/build.make src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o.provides.build
.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o.provides

src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o.provides.build: src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o


src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o: src/shared/CMakeFiles/amrl-shared-lib.dir/flags.make
src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o: ../src/shared/util/pthread_utils.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skm2019/Documents/CS393R/cs393r_starter/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o -c /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/pthread_utils.cc

src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.i"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/pthread_utils.cc > CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.i

src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.s"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/pthread_utils.cc -o CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.s

src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o.requires:

.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o.requires

src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o.provides: src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o.requires
	$(MAKE) -f src/shared/CMakeFiles/amrl-shared-lib.dir/build.make src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o.provides.build
.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o.provides

src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o.provides.build: src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o


src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o: src/shared/CMakeFiles/amrl-shared-lib.dir/flags.make
src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o: ../src/shared/util/timer.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skm2019/Documents/CS393R/cs393r_starter/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o -c /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/timer.cc

src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/amrl-shared-lib.dir/util/timer.cc.i"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/timer.cc > CMakeFiles/amrl-shared-lib.dir/util/timer.cc.i

src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/amrl-shared-lib.dir/util/timer.cc.s"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/timer.cc -o CMakeFiles/amrl-shared-lib.dir/util/timer.cc.s

src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o.requires:

.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o.requires

src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o.provides: src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o.requires
	$(MAKE) -f src/shared/CMakeFiles/amrl-shared-lib.dir/build.make src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o.provides.build
.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o.provides

src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o.provides.build: src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o


src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.o: src/shared/CMakeFiles/amrl-shared-lib.dir/flags.make
src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.o: ../src/shared/util/random.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skm2019/Documents/CS393R/cs393r_starter/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.o"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/amrl-shared-lib.dir/util/random.cc.o -c /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/random.cc

src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/amrl-shared-lib.dir/util/random.cc.i"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/random.cc > CMakeFiles/amrl-shared-lib.dir/util/random.cc.i

src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/amrl-shared-lib.dir/util/random.cc.s"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/random.cc -o CMakeFiles/amrl-shared-lib.dir/util/random.cc.s

src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.o.requires:

.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.o.requires

src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.o.provides: src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.o.requires
	$(MAKE) -f src/shared/CMakeFiles/amrl-shared-lib.dir/build.make src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.o.provides.build
.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.o.provides

src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.o.provides.build: src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.o


src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o: src/shared/CMakeFiles/amrl-shared-lib.dir/flags.make
src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o: ../src/shared/util/terminal_colors.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skm2019/Documents/CS393R/cs393r_starter/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o -c /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/terminal_colors.cc

src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.i"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/terminal_colors.cc > CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.i

src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.s"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skm2019/Documents/CS393R/cs393r_starter/src/shared/util/terminal_colors.cc -o CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.s

src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o.requires:

.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o.requires

src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o.provides: src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o.requires
	$(MAKE) -f src/shared/CMakeFiles/amrl-shared-lib.dir/build.make src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o.provides.build
.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o.provides

src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o.provides.build: src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o


# Object files for target amrl-shared-lib
amrl__shared__lib_OBJECTS = \
"CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o" \
"CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o" \
"CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o" \
"CMakeFiles/amrl-shared-lib.dir/util/random.cc.o" \
"CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o"

# External object files for target amrl-shared-lib
amrl__shared__lib_EXTERNAL_OBJECTS =

../lib/libamrl-shared-lib.so: src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o
../lib/libamrl-shared-lib.so: src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o
../lib/libamrl-shared-lib.so: src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o
../lib/libamrl-shared-lib.so: src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.o
../lib/libamrl-shared-lib.so: src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o
../lib/libamrl-shared-lib.so: src/shared/CMakeFiles/amrl-shared-lib.dir/build.make
../lib/libamrl-shared-lib.so: src/shared/CMakeFiles/amrl-shared-lib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/skm2019/Documents/CS393R/cs393r_starter/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX shared library ../../../lib/libamrl-shared-lib.so"
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/amrl-shared-lib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/shared/CMakeFiles/amrl-shared-lib.dir/build: ../lib/libamrl-shared-lib.so

.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/build

src/shared/CMakeFiles/amrl-shared-lib.dir/requires: src/shared/CMakeFiles/amrl-shared-lib.dir/util/helpers.cc.o.requires
src/shared/CMakeFiles/amrl-shared-lib.dir/requires: src/shared/CMakeFiles/amrl-shared-lib.dir/util/pthread_utils.cc.o.requires
src/shared/CMakeFiles/amrl-shared-lib.dir/requires: src/shared/CMakeFiles/amrl-shared-lib.dir/util/timer.cc.o.requires
src/shared/CMakeFiles/amrl-shared-lib.dir/requires: src/shared/CMakeFiles/amrl-shared-lib.dir/util/random.cc.o.requires
src/shared/CMakeFiles/amrl-shared-lib.dir/requires: src/shared/CMakeFiles/amrl-shared-lib.dir/util/terminal_colors.cc.o.requires

.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/requires

src/shared/CMakeFiles/amrl-shared-lib.dir/clean:
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared && $(CMAKE_COMMAND) -P CMakeFiles/amrl-shared-lib.dir/cmake_clean.cmake
.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/clean

src/shared/CMakeFiles/amrl-shared-lib.dir/depend:
	cd /home/skm2019/Documents/CS393R/cs393r_starter/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/skm2019/Documents/CS393R/cs393r_starter /home/skm2019/Documents/CS393R/cs393r_starter/src/shared /home/skm2019/Documents/CS393R/cs393r_starter/build /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared /home/skm2019/Documents/CS393R/cs393r_starter/build/src/shared/CMakeFiles/amrl-shared-lib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/shared/CMakeFiles/amrl-shared-lib.dir/depend

