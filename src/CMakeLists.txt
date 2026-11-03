set(BINARY ${CMAKE_PROJECT_NAME})

add_library(${BINARY}_lib STATIC
        Dist.h Dist.cpp
        DistributionFactory.h DistributionFactory.cpp
        Exceptions.h macros.h troolean.h troolean.cpp
        Assumptions.h Assumptions.cpp)
target_link_libraries(${BINARY}_lib PUBLIC ${GINAC_LIBRARIES})

if (PRODIGY_GINAC_BUILD_DEMO)
    add_executable(${BINARY}_demo demo.cpp)
    target_link_libraries(${BINARY}_demo PUBLIC ${BINARY}_lib)
endif()
