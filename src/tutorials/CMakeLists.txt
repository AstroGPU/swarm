CUDA_ADD_EXECUTABLE(simple simple.cpp ../utils/utils.cpp)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(simple swarm ${Boost_LIBRARIES})

