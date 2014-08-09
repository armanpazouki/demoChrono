#--------------------------------------------------------------
#
# Example of file to be used in CMake to build an external 
# project based on Chrono::Engine
#
#--------------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)


#--------------------------------------------------------------
#         Modify the project name if you want: 

PROJECT(demo_bricks)



#--------------------------------------------------------------
# NOTE!   use find_package() to define ChronoEngine libraries and variables, 
#         along with some optional units (components).
#         This requires that you copied FindChronoEngine.cmake in
#         the Modules/ directory of your Cmake installation. If it 
#         is not there, you can copy it in a cmake/ directory in your project, 
#         and add that directory to the search path of CMake by typing
#         set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(ChronoEngine   COMPONENTS unit_IRRLICHT )



#--------------------------------------------------------------
#         After the ChronoEngine package has been found, you
#         can add its include directories with the headers. 
#         Use ${CHRONOENGINE_INCLUDES} for this. 

INCLUDE_DIRECTORIES(  ${CHRONOENGINE_INCLUDES}  )
#INCLUDE_DIRECTORIES(  
#${CH_SDKDIR}
#${CH_SDKDIR}/collision/bullet
#"/usr/include/irrlicht"
#)

#--------------------------------------------------------------
#         Set c++ sources for building the exe, as usual in CMake

# Add executables


SET(DEMOS ${DEMOS}
    demo_bricks
    demo_ice1
    demo_ice2
    demo_excavator
)

FOREACH(PROGRAM ${DEMOS})
	MESSAGE(STATUS "...add ${PROGRAM}")

	ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
	SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

	#--------------------------------------------------------------
	#         The following is needed for some platforms
	SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES 
	LINK_FLAGS "${CH_LINKERFLAG_EXE}" )

	#--------------------------------------------------------------
	#         Uncomment the following if you want the executable to 
	#         go into some specific directory: 
	#            SET (CH_BINDIR_DEBUG   "" CACHE PATH "Where are your Chrono debug binaries (.dll, .exe etc.) installed?")
	#            SET (CH_BINDIR_RELEASE "" CACHE PATH "Where are your Chrono release binaries (.dll, .exe etc.) installed?")
	#            set_target_property(demo_bricks RUNTIME_OUTPUT_DIRECTORY_DEBUG   "${CH_BINDIR_DEBUG}")
	#            set_target_property(demo_bricks RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CH_BINDIR_RELEASE}")
	#         This done, you do not need to copy all the ChronoEngine.dll 
	#         files in the output directory; also, the bin/data/ dir
	#         will be shared by all programs.
	# OTHERWISE, VICEVERSA, copy the .dll files in your .exe folder.


	#--------------------------------------------------------------
	#  NOTE!  In the following, at least add the following line
	#         in order to link all the required libraries (the main 
	#         ChronoEngine library and all other libraries of units 
	#         defined with the COMPONENTS parameters in find_package)
	#          

	TARGET_LINK_LIBRARIES(${PROGRAM}
	${CHRONOENGINE_LIBRARIES}
	)









#    ADD_DEPENDENCIES(${PROGRAM} ${LIBRARIES})

 #   INSTALL(TARGETS ${PROGRAM} DESTINATION bin)

ENDFOREACH(PROGRAM)

















# this is only to simplify the GUI of CMake 
MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)
