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
CMAKE_SOURCE_DIR = /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/build

# Include any dependencies generated for this target.
include CMakeFiles/DSO_Course_1.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/DSO_Course_1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/DSO_Course_1.dir/flags.make

CMakeFiles/DSO_Course_1.dir/src/main.cpp.o: CMakeFiles/DSO_Course_1.dir/flags.make
CMakeFiles/DSO_Course_1.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/DSO_Course_1.dir/src/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DSO_Course_1.dir/src/main.cpp.o -c /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/src/main.cpp

CMakeFiles/DSO_Course_1.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DSO_Course_1.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/src/main.cpp > CMakeFiles/DSO_Course_1.dir/src/main.cpp.i

CMakeFiles/DSO_Course_1.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DSO_Course_1.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/src/main.cpp -o CMakeFiles/DSO_Course_1.dir/src/main.cpp.s

CMakeFiles/DSO_Course_1.dir/src/main.cpp.o.requires:

.PHONY : CMakeFiles/DSO_Course_1.dir/src/main.cpp.o.requires

CMakeFiles/DSO_Course_1.dir/src/main.cpp.o.provides: CMakeFiles/DSO_Course_1.dir/src/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/DSO_Course_1.dir/build.make CMakeFiles/DSO_Course_1.dir/src/main.cpp.o.provides.build
.PHONY : CMakeFiles/DSO_Course_1.dir/src/main.cpp.o.provides

CMakeFiles/DSO_Course_1.dir/src/main.cpp.o.provides.build: CMakeFiles/DSO_Course_1.dir/src/main.cpp.o


CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o: CMakeFiles/DSO_Course_1.dir/flags.make
CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o: ../src/geometry.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o -c /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/src/geometry.cpp

CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/src/geometry.cpp > CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.i

CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/src/geometry.cpp -o CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.s

CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o.requires:

.PHONY : CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o.requires

CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o.provides: CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o.requires
	$(MAKE) -f CMakeFiles/DSO_Course_1.dir/build.make CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o.provides.build
.PHONY : CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o.provides

CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o.provides.build: CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o


CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o: CMakeFiles/DSO_Course_1.dir/flags.make
CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o: ../src/photometric.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o -c /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/src/photometric.cpp

CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/src/photometric.cpp > CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.i

CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/src/photometric.cpp -o CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.s

CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o.requires:

.PHONY : CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o.requires

CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o.provides: CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o.requires
	$(MAKE) -f CMakeFiles/DSO_Course_1.dir/build.make CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o.provides.build
.PHONY : CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o.provides

CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o.provides.build: CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o


# Object files for target DSO_Course_1
DSO_Course_1_OBJECTS = \
"CMakeFiles/DSO_Course_1.dir/src/main.cpp.o" \
"CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o" \
"CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o"

# External object files for target DSO_Course_1
DSO_Course_1_EXTERNAL_OBJECTS =

../bin/DSO_Course_1: CMakeFiles/DSO_Course_1.dir/src/main.cpp.o
../bin/DSO_Course_1: CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o
../bin/DSO_Course_1: CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o
../bin/DSO_Course_1: CMakeFiles/DSO_Course_1.dir/build.make
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_shape.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_stitching.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_superres.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_videostab.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_aruco.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_bgsegm.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_bioinspired.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_ccalib.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_datasets.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_dpm.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_face.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_freetype.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_fuzzy.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_hdf.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_line_descriptor.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_optflow.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_plot.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_reg.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_saliency.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_stereo.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_structured_light.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_surface_matching.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_text.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_ximgproc.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_xobjdetect.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_xphoto.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_video.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_viz.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_phase_unwrapping.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_rgbd.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_calib3d.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_features2d.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_flann.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_objdetect.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_ml.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_highgui.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_photo.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_videoio.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_imgcodecs.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_imgproc.so.3.2.0
../bin/DSO_Course_1: /usr/lib/x86_64-linux-gnu/libopencv_core.so.3.2.0
../bin/DSO_Course_1: CMakeFiles/DSO_Course_1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable ../bin/DSO_Course_1"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DSO_Course_1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/DSO_Course_1.dir/build: ../bin/DSO_Course_1

.PHONY : CMakeFiles/DSO_Course_1.dir/build

CMakeFiles/DSO_Course_1.dir/requires: CMakeFiles/DSO_Course_1.dir/src/main.cpp.o.requires
CMakeFiles/DSO_Course_1.dir/requires: CMakeFiles/DSO_Course_1.dir/src/geometry.cpp.o.requires
CMakeFiles/DSO_Course_1.dir/requires: CMakeFiles/DSO_Course_1.dir/src/photometric.cpp.o.requires

.PHONY : CMakeFiles/DSO_Course_1.dir/requires

CMakeFiles/DSO_Course_1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/DSO_Course_1.dir/cmake_clean.cmake
.PHONY : CMakeFiles/DSO_Course_1.dir/clean

CMakeFiles/DSO_Course_1.dir/depend:
	cd /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1 /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1 /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/build /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/build /home/jie/Desktop/shenlan/orbslam/slam_course/DSO/HW1/build/CMakeFiles/DSO_Course_1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/DSO_Course_1.dir/depend

