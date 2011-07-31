CUDA_ADD_EXECUTABLE(simple simple.cpp ../utils/utils.cpp)
TARGET_LINK_LIBRARIES(simple swarm ${Boost_LIBRARIES})

