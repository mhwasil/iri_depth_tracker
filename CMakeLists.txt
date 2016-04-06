cmake_minimum_required(VERSION 2.4.6)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries
set(ROS_BUILD_TYPE Release)

rosbuild_init()

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

#uncomment if you have defined messages
#rosbuild_genmsg()
#uncomment if you have defined services
#rosbuild_gensrv()

#common commands for building c++ executables and libraries
#rosbuild_add_library(${PROJECT_NAME} src/example.cpp)
#target_link_libraries(${PROJECT_NAME} another_library)
#rosbuild_add_boost_directories()
#rosbuild_link_boost(${PROJECT_NAME} thread)
#rosbuild_add_executable(example examples/example.cpp)
#target_link_libraries(example ${PROJECT_NAME})set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -std=c++0x -fopenmp -fPIC")

rosbuild_add_library(lmath src/chol.cpp src/compute_prob1.cpp
  src/crossTimesMatrix.cpp src/cumulative_sum.cpp src/diag.cpp src/eig.cpp
   src/estimateRigidTransform.cpp src/exp.cpp src/expm.cpp src/eye.cpp src/gradient.cpp
    src/image_warping.cpp src/init_variables.cpp src/log.cpp src/mean.cpp src/meshgrid.cpp src/mldivide.cpp src/mrdivide.cpp
      src/norm.cpp src/Optimal_affine_tracking_3d16_fast_realtime_data.cpp
       src/Optimal_affine_tracking_3d16_fast_realtime_initialize.cpp 
       src/Optimal_affine_tracking_3d16_fast_realtime_rtwutil.cpp 
       src/Optimal_affine_tracking_3d16_fast_realtime_terminate.cpp src/quat2rot.cpp 
       src/rand.cpp src/randn.cpp src/repmat.cpp src/resampling.cpp src/round.cpp src/rtGetNaN.cpp
        src/rtGetInf.cpp src/rt_nonfinite.cpp src/sqrt.cpp src/sum.cpp src/svd.cpp)



rosbuild_add_executable(depth_tracker src/main.cpp)
target_link_libraries(depth_tracker lmath)
rosbuild_add_openmp_flags(depth_tracker)
