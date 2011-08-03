CUDA_ADD_EXECUTABLE(stability_test stability_test.cpp utils.cpp)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(stability_test swarm ${Boost_LIBRARIES})

ADD_EXECUTABLE(test_energy test_energy.cpp)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(test_energy swarm)

ADD_EXECUTABLE(snapshot2text snapshot2text.cpp utils.cpp ../swarm/snapshot.cpp)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(snapshot2text ${Boost_LIBRARIES})

ADD_EXECUTABLE(text2snapshot text2snapshot.cpp utils.cpp ../swarm/snapshot.cpp)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(text2snapshot ${Boost_LIBRARIES})

