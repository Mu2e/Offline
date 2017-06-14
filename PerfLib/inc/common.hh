#ifndef _perflib_common_h_
#define _perflib_common_h_

#include <fcntl.h>
#include <unistd.h>
#include <chrono>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <list>
#include <memory>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <vector>
//#include <boost/exception/diagnostic_information.hpp>
#include <cstring>
namespace perflib {
struct result_t {
  int code;
  std::string message;
};
}

#endif /* _perflib_common_h_ */
