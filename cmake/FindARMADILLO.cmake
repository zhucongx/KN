set(FIND_ARMADILLO_PATHS ~/armadillo /usr/ /usr/local)

find_path(ARMADILLO_INCLUDE_DIRS
        NAMES armadillo
        PATH_SUFFIXES include
        PATHS ${FIND_ARMADILLO_PATHS})
find_library(ARMADILLO_LIBRARIES
        NAMES armadillo
        HINTS ${FIND_ARMADILLO_PATHS})
