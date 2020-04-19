if (APPLE)
    message("-- Looking for armadillo on macOS")
    set(FIND_ARMADILLO_PATHS /usr/ /usr/local)
elseif (UNIX)
    message("-- Looking for armadillo on Linux")
    set(FIND_ARMADILLO_PATHS /usr/ ~/armadillo/)
endif()

find_path(ARMADILLO_INCLUDE_DIRS
        NAMES armadillo
        PATH_SUFFIXES include
        PATHS ${FIND_ARMADILLO_PATHS})
find_library(ARMADILLO_LIBRARIES
        NAMES armadillo
        PATH_SUFFIXES lib64 lib
        PATHS ${FIND_ARMADILLO_PATHS})