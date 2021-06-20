cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Catch2 REQUIRED)
find_package(HTSlib REQUIRED)

include_directories(nlohmann_json::nlohmann_json)

enable_testing()

add_custom_target(bam_pileup_data COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/make-pileup.sh ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_gtf test_gtf.cpp)
target_link_libraries(test_gtf PRIVATE Catch2::Catch2)
add_test(NAME test_gtf COMMAND test_gtf WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(stabby_test_cons stabby/test_cons.cpp)
target_link_libraries(stabby_test_cons PRIVATE sdsl divsufsort divsufsort64 Catch2::Catch2)
add_test(NAME stabby_test_cons COMMAND stabby_test_cons)

add_executable(stabby_test_regressions stabby/test_regressions.cpp)
target_link_libraries(stabby_test_regressions PRIVATE sdsl divsufsort divsufsort64 nlohmann_json::nlohmann_json Catch2::Catch2)
add_test(NAME stabby_test_regressions COMMAND stabby_test_regressions WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

#add_executable(vcf_test_vcf_parser vcf/test_vcf_parser.cpp)
#target_link_libraries(vcf_test_vcf_parser PRIVATE Catch2::Catch2)
#add_test(NAME vcf_test_vcf_parser  COMMAND vcf_test_vcf_parser)

add_executable(vcf_test_vcf vcf/test_vcf.cpp)
target_link_libraries(vcf_test_vcf hts Catch2::Catch2 Threads::Threads)
add_test(NAME vcf_test_vcf  COMMAND vcf_test_vcf WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(bam_test_bam test_bam.cpp)
target_link_libraries(bam_test_bam hts Catch2::Catch2 Threads::Threads)
add_dependencies(bam_test_bam bam_pileup_data)
add_test(NAME bam_test_bam  COMMAND bam_test_bam WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

