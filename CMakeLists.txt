cmake_minimum_required(VERSION 3.14)                                                  
project(rtsim CXX)                                                                    
                                                                                          
set(CMAKE_CXX_STANDARD 17)                                                            
set(CMAKE_CXX_STANDARD_REQUIRED ON)                                                   

add_compile_options(-O3 -march=native)                                

add_library(rtsim_lib OBJECT src/rtsim.cc) # Define rtsim object library                                                                                          

add_executable(sim src/main.cc) # Create and compile main as sim target                
target_link_libraries(sim PRIVATE rtsim_lib) # Link sim to depend on rtsim                            

add_executable(reg_test src/regression_test.cc) # Create and compile regression_test code as reg_test target 
target_link_libraries(reg_test PRIVATE rtsim_lib) # Link reg_test target to depend on rtsim

enable_testing() # Activate CTest framework
add_test(NAME RegressionTest COMMAND reg_test) # Register reg_test target as a test