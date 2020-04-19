if (APPLE)
    message("-- Looking for armadillo on macOS")
    set(FIND_ARMADILLO_PATHS /usr/local/)
elseif (UNIX)
    message("-- Looking for armadillo on Linux")
    set(FIND_ARMADILLO_PATHS /usr/ ~/armadillo/)
endif()

find_path(ARMADILLO_INCLUDE_DIRS armadillo
        PATH_SUFFIXES include
        PATHS ${FIND_ARMADILLO_PATHS})
find_library(ARMADILLO_LIBRARIES
        NAMES libarmadillo
        PATH_SUFFIXES lib64
        PATHS ${FIND_ARMADILLO_PATHS})