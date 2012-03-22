FIND_PACKAGE(PythonLibs)
FIND_PACKAGE(Boost COMPONENTS PYTHON)
IF(Boost_PYTHON_FOUND)
	IF(PYTHONLIBS_FOUND)
		SWARM_ADD_LIBRARY(swarmng_ext SHARED
			module.cpp
		)
		INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
		TARGET_LINK_LIBRARIES(swarmng_ext ${Python_LIBRARIES})
	ELSE()
		MESSAGE("Python headers were not found")
	endif()
else()
	MESSAGE("Boost::Python is not found")
ENDIF()

