set(BINARY ${CMAKE_PROJECT_NAME})

add_executable(${BINARY}_demo demo.cpp)

add_library(${BINARY}_lib STATIC
        Formula.h Formula.cpp
        ParametricDist.h ParametricDist.cpp
        DistributionFactory.h DistributionFactory.cpp
        Exceptions.h)
target_link_libraries(${BINARY}_lib PUBLIC ${GINAC_LIBRARIES})

target_link_libraries(${BINARY}_demo PUBLIC ${BINARY}_lib)
