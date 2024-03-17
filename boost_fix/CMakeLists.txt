cet_make_library(INTERFACE INSTALLED_PATH_BASE Offline
    SOURCE
      accumulators/statistics.hpp
      accumulators/statistics/stats.hpp
    LIBRARIES INTERFACE
)

install_headers(USE_PROJECT_NAME SUBDIRS accumulators accumulators/statistics)
