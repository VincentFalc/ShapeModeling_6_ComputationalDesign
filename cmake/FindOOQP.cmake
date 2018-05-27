# TRY TO FIND THE INCLUDE DIRECTORY
find_path(OOQP_INCLUDE_DIR
  QpGenData.h
  HINTS 
    ${CMAKE_SOURCE_DIR}/../OOQP/include/ooqp
	$ENV{OOQP}/include/ooqp
    /usr/local/include/ooqp/ 
    ${PROJECT_SOURCE_DIR}/../../libs/thirdPartyCode/OOQP/include
)

# BLAS
find_library(BLAS_LIBRARIES
  NAMES 
    blas
  HINTS
    ${CMAKE_SOURCE_DIR}/../OOQP/lib
    /usr/lib/x86_64-linux-gnu/
    /usr/local/libs/ 
    ${CMAKE_SOURCE_DIR}/../libs/thirdPartyCode/CLAPACK/BLAS/Release
)
if(BLAS_LIBRARIES)
  message(STATUS "Found BLAS libraries:" ${BLAS_LIBRARIES})
else()
  message(FATAL_ERROR "Could not find BLAS libraries.")
endif()

# Dependent packages, BLAS and HSL
# set(OOQP_LIBRARIES)
# find_package(BLAS REQUIRED)
# if(BLAS_FOUND)
#   message("BLAS FOUND")
# else()
#   message(STATUS "OOQP requires BLAS")
# endif()

# TODO: do we need these libs on windows?
if (WIN32)
  # I77
  set(I77_LIB I77)
  if(WIN32)
    set(I77_LIB lib${I77_LIB})
  endif()
  find_library(F2CLIBS_I77
    NAMES
      ${I77_LIB}
    HINTS
	  ${CMAKE_SOURCE_DIR}/../OOQP/lib
      /usr/local/lib/
      /usr/lib/
      /usr/lib/x86_64-linux-gnu/
      ${CMAKE_SOURCE_DIR}/../libs/thirdPartyCode/CLAPACK/F2CLIBS/ReleaseI77
  )

  # I77
  set(F77_LIB F77)
  if(WIN32)
    set(F77_LIB lib${F77_LIB})
  endif()
  find_library(F2CLIBS_F77
    NAMES
      ${F77_LIB}
    HINTS
	  ${CMAKE_SOURCE_DIR}/../OOQP/lib
      /usr/local/lib/
      ${CMAKE_SOURCE_DIR}/../libs/thirdPartyCode/CLAPACK/F2CLIBS/ReleaseF77
  )
  set(F2CLIBS_LIBRARIES ${F2CLIBS_I77} ${F2CLIBS_F77})
  if(F2CLIBS_LIBRARIES)
    message(STATUS "Found F2CLIBS libraries:" ${F2CLIBS_LIBRARIES})
  else()
    message(FATAL_ERROR "Could not find F2CLIBS libraries.:" ${F2CLIBS_LIBRARIES})
  endif()
endif(WIN32)

# CLAPACK
set(CLAPCK_LIB lapack)
if(WIN32)
  set(CLAPCK_LIB c${CLAPCK_LIB})
endif()
find_library(CLAPACK_LIBRARIES
  NAMES 
    ${CLAPCK_LIB}
  HINTS 
    ${CMAKE_SOURCE_DIR}/../OOQP/lib
    /usr/lib/x86_64-linux-gnu/
    /usr/local/lib/
    /usr/lib/
    ${CMAKE_SOURCE_DIR}/../libs/thirdPartyCode/CLAPACK/Release
)
if(CLAPACK_LIBRARIES)
  message(STATUS "Found CLAPACK libraries:" ${CLAPACK_LIBRARIES})
else()
  message(FATAL_ERROR "Could not find CLAPACK libraries.")
endif()

set(OOQP_LIBRARIES ${BLAS_LIBRARIES} ${F2CLIBS_LIBRARIES} ${CLAPACK_LIBRARIES})

# find_package(BLAS REQUIRED)
################################

# Dependent packages, BLAS and HSL
# set(OOQP_LIBRARIES)
# find_package(BLAS REQUIRED)
# if(BLAS_FOUND)
#   message("BLAS FOUND")
#   set(OOQP_LIBRARIES ${OOQP_LIBRARIES} ${BLAS_LIBRARIES})
# else()
#   message(STATUS "OOQP requires BLAS")
# endif()
# find_package(HSL QUIET)
# if(HSL_FOUND)
#   set(OOQP_LIBRARIES ${OOQP_LIBRARIES} ${HSL_LIBRARIES})
# else()
#   message(STATUS "OOQP requires HSL")
# endif()

if(OOQP_INCLUDE_DIR)
  set(OOQP_FOUND_INCLUDE TRUE)
  set(OOQP_INCLUDE_DIRS
  ${OOQP_INCLUDE_DIR})
  message(STATUS "Found OOQP include dirs: ${OOQP_INCLUDE_DIRS}")
else()
  message(STATUS "Could not find OOQP include dir")
endif()

# TRY TO FIND THE LIBRARIES
set(OOQP_LIBS_LIST
  ooqpgensparse ooqpsparse ooqpgondzio ooqpbase ma27
)

set(OOQP_FOUND_LIBS TRUE)
foreach(LIB ${OOQP_LIBS_LIST})
  
  if(WIN32)
    set(LIB lib${LIB}.lib)
  endif()

  find_library(OOQP_LIB_${LIB}
    NAMES ${LIB}
    HINTS 
	  ${CMAKE_SOURCE_DIR}/../OOQP/lib
      /usr/local/libs/ 
      ${PROJECT_SOURCE_DIR}/../../libs/thirdPartyCode/OOQP/lib/
  )
  if(OOQP_LIB_${LIB})
    set(OOQP_LIBRARIES ${OOQP_LIBRARIES} ${OOQP_LIB_${LIB}})
  else()
    message(FATAL_ERROR "Could not find " ${LIB})
    set(OOQP_FOUND_LIBS FALSE)
  endif()
endforeach()

# TODO: this is not very clean, use find package
if (UNIX)
  set(OOQP_LIBRARIES ${OOQP_LIBRARIES} gfortran blas)
endif (UNIX)


# print OOQP_LIBRARIES
if(OOQP_FOUND_LIBS)
  message(STATUS "Found OOQP libraries: ${OOQP_LIBRARIES}")
else()
  message(STATUS "Cound not find OOQP libraries")
endif()

# SUCCESS if BOTH THE LIBRARIES AND THE INCLUDE DIRECTORIES WERE FOUND
if(OOQP_FOUND_INCLUDE AND OOQP_FOUND_LIBS AND BLAS_FOUND AND HSL_FOUND)
  set(OOQP_FOUND TRUE)
  message(STATUS "Found OOQP")
elseif()
  message(STATUS "Cound not find OOQP")
endif()
