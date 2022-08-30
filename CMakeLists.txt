cmake_minimum_required(VERSION 3.2)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#add_compile_options(-msse4.2)

project(Scindo VERSION 1.0)
set(Scindo_VERSION 1.0)

include(CTest)

# Function to prepend the subdirectory to source files in subdirectories
FUNCTION(PREPEND var )
   SET(listVar "")
   FOREACH(f ${${var}})
       LIST(APPEND listVar "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

# After a normal build, we can specify the location of various outputs of the
# build. We put executables and static libraries outside the build directory in
# bin/ and lib/, respectively.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE RelWithDebInfo)

MESSAGE(DEBUG "${CMAKE_MODULE_PATH}")

find_package(Threads REQUIRED)
find_package(HTSlib REQUIRED)
find_package(ZLIB REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(docopt COMPONENTS CXX REQUIRED)
find_package(pegtl)
 
#SET(Boost_USE_STATIC_LIBS ON)
FIND_PACKAGE(Boost COMPONENTS log log_setup iostreams REQUIRED)

add_subdirectory(src)

add_subdirectory(include)

include_directories(${DOCOPT_INCLUDE_DIRS})
include_directories(include)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

add_executable(scindo ${Scindo_SRC})
target_link_libraries(
    scindo
    docopt
    hts
    sdsl divsufsort divsufsort64
    ${ZLIB_LIBRARIES}
    nlohmann_json::nlohmann_json
    Boost::log_setup
    ${Boost_LIBRARIES}
    Threads::Threads)

add_executable(hecil ${Hecil_SRC})
target_link_libraries(
    hecil
    docopt
    hts
    sdsl divsufsort divsufsort64
    ${ZLIB_LIBRARIES}
    nlohmann_json::nlohmann_json
    Boost::log_setup
    ${Boost_LIBRARIES}
    Threads::Threads)

add_executable(asya ${Asya_SRC})
target_link_libraries(
    asya
    docopt
    hts
    sdsl divsufsort divsufsort64
    ${ZLIB_LIBRARIES}
    nlohmann_json::nlohmann_json
    Boost::log_setup
    ${Boost_LIBRARIES}
    Threads::Threads)

add_executable(perya ${Perya_SRC})
target_link_libraries(
    perya
    docopt
    hts
    sdsl divsufsort divsufsort64
    ${ZLIB_LIBRARIES}
    nlohmann_json::nlohmann_json
    Boost::log_setup
    ${Boost_LIBRARIES}
    Threads::Threads)

add_executable(crispin ${Crispin_SRC})
target_link_libraries(
    crispin
    PRIVATE
    docopt
    hts
    sdsl divsufsort divsufsort64
    ${ZLIB_LIBRARIES}
    nlohmann_json::nlohmann_json
    Boost::log_setup
    ${Boost_LIBRARIES}
    Threads::Threads)

add_executable(fraggle ${Fraggle_SRC})
target_link_libraries(
    fraggle
    PRIVATE
    docopt
    sdsl divsufsort divsufsort64
    ${ZLIB_LIBRARIES}
    nlohmann_json::nlohmann_json
    Boost::log_setup
    ${Boost_LIBRARIES}
    Threads::Threads)
    

install(TARGETS scindo hecil perya crispin
        CONFIGURATIONS Release)
install(TARGETS scindo hecil perya crispin
        CONFIGURATIONS RelWithDebInfo)
