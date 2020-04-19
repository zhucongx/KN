if (APPLE)
    message("-- Looking for armadillo on macOS")
    set(FIND_ARMADILLO_PATHS "/usr/ /usr/local")
elseif (UNIX)
    message("-- Looking for armadillo on Linux")
    set(FIND_ARMADILLO_PATHS "/usr/ ~/armadillo/")
endif()

find_path(ARMADILLO_INCLUDE_DIRS
        NAMES armadillo
        PATH_SUFFIXES include
        PATHS ${FIND_ARMADILLO_PATHS})
find_library(ARMADILLO_LIBRARIES
        NAMES libarmadillo.dylib libarmadillo.so
        PATH_SUFFIXES lib lib64
        PATHS ${FIND_ARMADILLO_PATHS})