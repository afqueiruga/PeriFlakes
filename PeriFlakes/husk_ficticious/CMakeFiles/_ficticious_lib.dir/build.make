# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious

# Include any dependencies generated for this target.
include CMakeFiles/_ficticious_lib.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/_ficticious_lib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/_ficticious_lib.dir/flags.make

ficticious_swigPYTHON_wrap.c: ficticious_swig.i
ficticious_swigPYTHON_wrap.c: ficticious_swig.i
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Swig source"
	/opt/local/bin/cmake -E make_directory /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious
	/opt/local/bin/swig -python -outdir /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious -I/Users/afq/Documents/Research/opt/lblsuper2/include -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious -o /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/ficticious_swigPYTHON_wrap.c /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/ficticious_swig.i

ficticious_lib.py: ficticious_swigPYTHON_wrap.c
	@$(CMAKE_COMMAND) -E touch_nocreate ficticious_lib.py

CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o: CMakeFiles/_ficticious_lib.dir/flags.make
CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o: ficticious_swigPYTHON_wrap.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o"
	/opt/local/bin/gcc-mp-5 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o   -c /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/ficticious_swigPYTHON_wrap.c

CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.i"
	/opt/local/bin/gcc-mp-5 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/ficticious_swigPYTHON_wrap.c > CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.i

CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.s"
	/opt/local/bin/gcc-mp-5 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/ficticious_swigPYTHON_wrap.c -o CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.s

CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o.requires:

.PHONY : CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o.requires

CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o.provides: CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o.requires
	$(MAKE) -f CMakeFiles/_ficticious_lib.dir/build.make CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o.provides.build
.PHONY : CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o.provides

CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o.provides.build: CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o


CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o: CMakeFiles/_ficticious_lib.dir/flags.make
CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o: bobaru_n.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o"
	/opt/local/bin/gcc-mp-5 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS)  -fPIC -o CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o   -c /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/bobaru_n.c

CMakeFiles/_ficticious_lib.dir/bobaru_n.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/_ficticious_lib.dir/bobaru_n.c.i"
	/opt/local/bin/gcc-mp-5 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS)  -fPIC -E /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/bobaru_n.c > CMakeFiles/_ficticious_lib.dir/bobaru_n.c.i

CMakeFiles/_ficticious_lib.dir/bobaru_n.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/_ficticious_lib.dir/bobaru_n.c.s"
	/opt/local/bin/gcc-mp-5 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS)  -fPIC -S /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/bobaru_n.c -o CMakeFiles/_ficticious_lib.dir/bobaru_n.c.s

CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o.requires:

.PHONY : CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o.requires

CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o.provides: CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o.requires
	$(MAKE) -f CMakeFiles/_ficticious_lib.dir/build.make CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o.provides.build
.PHONY : CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o.provides

CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o.provides.build: CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o


CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o: CMakeFiles/_ficticious_lib.dir/flags.make
CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o: bobaru_y.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o"
	/opt/local/bin/gcc-mp-5 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS)  -fPIC -o CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o   -c /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/bobaru_y.c

CMakeFiles/_ficticious_lib.dir/bobaru_y.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/_ficticious_lib.dir/bobaru_y.c.i"
	/opt/local/bin/gcc-mp-5 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS)  -fPIC -E /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/bobaru_y.c > CMakeFiles/_ficticious_lib.dir/bobaru_y.c.i

CMakeFiles/_ficticious_lib.dir/bobaru_y.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/_ficticious_lib.dir/bobaru_y.c.s"
	/opt/local/bin/gcc-mp-5 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS)  -fPIC -S /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/bobaru_y.c -o CMakeFiles/_ficticious_lib.dir/bobaru_y.c.s

CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o.requires:

.PHONY : CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o.requires

CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o.provides: CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o.requires
	$(MAKE) -f CMakeFiles/_ficticious_lib.dir/build.make CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o.provides.build
.PHONY : CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o.provides

CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o.provides.build: CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o


CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o: CMakeFiles/_ficticious_lib.dir/flags.make
CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o: bobaru_n3.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o"
	/opt/local/bin/gcc-mp-5 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS)  -fPIC -o CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o   -c /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/bobaru_n3.c

CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.i"
	/opt/local/bin/gcc-mp-5 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS)  -fPIC -E /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/bobaru_n3.c > CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.i

CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.s"
	/opt/local/bin/gcc-mp-5 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS)  -fPIC -S /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/bobaru_n3.c -o CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.s

CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o.requires:

.PHONY : CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o.requires

CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o.provides: CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o.requires
	$(MAKE) -f CMakeFiles/_ficticious_lib.dir/build.make CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o.provides.build
.PHONY : CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o.provides

CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o.provides.build: CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o


# Object files for target _ficticious_lib
_ficticious_lib_OBJECTS = \
"CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o" \
"CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o" \
"CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o" \
"CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o"

# External object files for target _ficticious_lib
_ficticious_lib_EXTERNAL_OBJECTS =

_ficticious_lib.so: CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o
_ficticious_lib.so: CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o
_ficticious_lib.so: CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o
_ficticious_lib.so: CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o
_ficticious_lib.so: CMakeFiles/_ficticious_lib.dir/build.make
_ficticious_lib.so: /opt/local/lib/libgsl.dylib
_ficticious_lib.so: /opt/local/lib/libgslcblas.dylib
_ficticious_lib.so: /opt/local/lib/libpython2.7.dylib
_ficticious_lib.so: /Users/afq/Documents/Research/opt/lblsuper2/lib/libcornflakes.a
_ficticious_lib.so: CMakeFiles/_ficticious_lib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking C shared module _ficticious_lib.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/_ficticious_lib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/_ficticious_lib.dir/build: _ficticious_lib.so

.PHONY : CMakeFiles/_ficticious_lib.dir/build

CMakeFiles/_ficticious_lib.dir/requires: CMakeFiles/_ficticious_lib.dir/ficticious_swigPYTHON_wrap.c.o.requires
CMakeFiles/_ficticious_lib.dir/requires: CMakeFiles/_ficticious_lib.dir/bobaru_n.c.o.requires
CMakeFiles/_ficticious_lib.dir/requires: CMakeFiles/_ficticious_lib.dir/bobaru_y.c.o.requires
CMakeFiles/_ficticious_lib.dir/requires: CMakeFiles/_ficticious_lib.dir/bobaru_n3.c.o.requires

.PHONY : CMakeFiles/_ficticious_lib.dir/requires

CMakeFiles/_ficticious_lib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/_ficticious_lib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/_ficticious_lib.dir/clean

CMakeFiles/_ficticious_lib.dir/depend: ficticious_swigPYTHON_wrap.c
CMakeFiles/_ficticious_lib.dir/depend: ficticious_lib.py
	cd /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious /Users/afq/Documents/Dropbox/LBL/PeriFlakes/PeriFlakes/husk_ficticious/CMakeFiles/_ficticious_lib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/_ficticious_lib.dir/depend

