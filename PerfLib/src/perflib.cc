#include "PerfLib/inc/perflib.hh"

// include <linux/hw_breakpoint.h>
#include <linux/perf_event.h>
#include <pthread.h>
#include <sched.h>
#include <sqlite3.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <utmpx.h>

using perf::PerfStats;

int PerfStats::fuse_core(int core_id) {
  if (core_id == -1) core_id = sched_getcpu();

  int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
  if (core_id < 0 || core_id >= num_cores) return EINVAL;

  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(core_id, &cpuset);

  pthread_t current_thread = pthread_self();
  int retval = pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);

  return retval == 0 ? core_id : -1;
}

void PerfStats::start_hwctr(size_t idx) {
  if (_conters_dec[idx].method == RDTSC) return;

  auto& desc = _conters_dec[idx];

  memset(&desc.attr, 0, sizeof(decltype(desc.attr)));
  desc.attr.type = desc.type;
  desc.attr.size = sizeof(decltype(desc.attr));
  desc.attr.config = desc.event_mask;
  desc.attr.inherit = 1;
  desc.attr.disabled = 1;
  desc.attr.exclude_kernel = 1;
  desc.attr.exclude_hv = 1;

  desc.file = syscall(__NR_perf_event_open, &desc.attr, 0, -1, -1, 0);

  if (desc.file < 0) std::cout << "Failed to start " << desc.name << ".\n";

  if (ioctl(desc.file, PERF_EVENT_IOC_RESET, 0) == -1) std::cout << "Failed to reset counter " << desc.name << ".\n";

  if (ioctl(desc.file, PERF_EVENT_IOC_ENABLE, 0) == -1) std::cout << "Failed to enable counter " << desc.name << ".\n";
}

void PerfStats::stop_hwctr(size_t idx) {
  if (_conters_dec[idx].method == RDTSC) return;

  auto& desc = _conters_dec[idx];

  if (ioctl(desc.file, PERF_EVENT_IOC_DISABLE, 0) == -1)
    std::cout << "Failed to disable counter " << desc.name << ".\n";

  close(desc.file);
  desc.file = -1;
}

void PerfStats::add_counter(unsigned int type, unsigned int event_mask, int method, std::string const& name,
                            unsigned int rdpmc_eax) {
  if (_counter_count == MAX_COUNTERS) throw std::runtime_error("Error: Counter array is full.");

  auto& desc = _conters_dec[_counter_count];
  desc.type = type;
  desc.event_mask = event_mask;
  desc.name = name;
  desc.method = method;
  desc.rdpmc_eax = rdpmc_eax;

  ++_counter_count;
}

void PerfStats::record_counters() {
  if (_sample_count < 2) return;

  sqlite3* db;
  std::stringstream sql;

  if (sqlite3_open("perfstats.sqlite3", &db) != SQLITE_OK) {
    std::cout << "Can't open database \n";
    goto close;
  }

  sql << "CREATE TABLE if not exists perfstats(\
  id INTEGER PRIMARY KEY, \
  timestamp DATETIME DEFAULT (STRFTIME(\'%Y-%m-%d %H:%M:%f\', \'now\', \'localtime\')),\
  observation TEXT, \
  iteration INTEGER,  \
  name TEXT, \
  value INTEGER \
  );\n";

  sql << "INSERT INTO perfstats ( observation, iteration, name, value) VALUES";

  for (size_t ismpl = 0; ismpl < _sample_count; ++ismpl, ++_written_sample_count) {
    for (size_t idx = 0; idx < _counter_count; ++idx) {
      sql << "(";
      sql << "\'" << _observation_name << "\'" << ',';
      sql << _written_sample_count << ',';
      sql << "\'" << _conters_dec[idx].name << "\'" << ',';
      sql << _samples[ismpl][idx];
      sql << "),\n";
    }
  }

  _sample_count = 0;

  sql.seekp(-2, sql.cur);
  sql << ";\n";

  if (sqlite3_exec(db, sql.str().c_str(), 0, 0, 0) != SQLITE_OK) std::cout << "SQL error:  sql=" << sql.str() << "\n";

close:
  sqlite3_close(db);
}
/*
void PerfStats::record_counters(std::string const& observation, size_t iteration) {
  if(iteration==0)
    return;

  sqlite3* db;
  std::stringstream sql;

  if (sqlite3_open("perfstats.sqlite3", &db) != SQLITE_OK) {
    std::cout << "Can't open database \n";
    goto close;
  }

  sql << "CREATE TABLE if not exists perfstats(\
  id INTEGER PRIMARY KEY, \
  timestamp DATETIME DEFAULT (STRFTIME(\'%Y-%m-%d %H:%M:%f\', \'now\', \'localtime\')),\
  observation TEXT, \
  iteration INTEGER,  \
  name TEXT, \
  value INTEGER \
  );\n";

  sql << "INSERT INTO perfstats ( observation, iteration, name, value) VALUES";

  for (size_t idx = 0; idx < _counter_count; ++idx) {
    sql << "(";
    sql << "\'" << observation << "\'" << ',';
    sql << iteration << ',';
    sql << "\'" << _conters_dec[idx].name << "\'" << ',';
    sql << _sample_count[iteration][idx];
    sql << "),\n";
  }
  sql.seekp(-2, sql.cur);
  sql << ";\n";

  if (sqlite3_exec(db, sql.str().c_str(), 0, 0, 0) != SQLITE_OK) std::cout << "SQL error:  sql=" << sql.str() << "\n";

close:
  sqlite3_close(db);
}

*/

void PerfStats::record_metrics(metrics_t const& metrics, std::string const& observation, unsigned int iteration) {
  if (iteration == 0) return;

  sqlite3* db;
  std::stringstream sql;

  if (sqlite3_open("perfstats.sqlite3", &db) != SQLITE_OK) {
    std::cout << "Can't open database \n";
    goto close;
  }

  sql << "CREATE TABLE if not exists metrics(\
  id INTEGER PRIMARY KEY, \
  timestamp DATETIME DEFAULT (STRFTIME(\'%Y-%m-%d %H:%M:%f\', \'now\', \'localtime\')),\
  observation TEXT, \
  iteration INTEGER,  \
  name TEXT, \
  value REAL \
  );\n";

  sql << "INSERT INTO metrics ( observation, iteration, name, value) VALUES";

  for (auto const& metric : metrics) {
    sql << "(";
    sql << "\'" << observation << "\'" << ',';
    sql << iteration << ',';
    sql << "\'" << metric.first << "\'" << ',';
    sql << metric.second;
    sql << "),\n";
  }

  sql.seekp(-2, sql.cur);
  sql << ";\n";

  if (sqlite3_exec(db, sql.str().c_str(), 0, 0, 0) != SQLITE_OK) std::cout << "SQL error:  sql=" << sql.str() << "\n";

close:
  sqlite3_close(db);
}

static int report_counters_callback(void* NotUsed, int argc, char** argv, char** azColName) {
  std::stringstream* result = (std::stringstream*)NotUsed;

  for (int x = 0; x < argc; x += 2) *result << "  " << argv[x] << " = " << argv[x + 1] << "\n";

  return SQLITE_OK;
}

std::string PerfStats::report_counters() {
  sqlite3* db;
  std::stringstream sql, result;
  result << "Counters:\n";

  if (sqlite3_open("perfstats.sqlite3", &db) != SQLITE_OK) {
    std::cout << "Can't open database \n";
    goto close;
  }

  sql << "select name, value from (select * from perfstats order by id DESC limit " << _counter_count << ");";

  if (sqlite3_exec(db, sql.str().c_str(), report_counters_callback, (void*)&result, 0) != SQLITE_OK)
    std::cout << "SQL error:  sql=" << sql.str() << "\n";

close:
  result << "\n";
  sqlite3_close(db);

  return result.str();
}

void PerfStats::init() {
  _begin.fill(0);
  _end.fill(0);
  _delta.fill(0);

  auto core_id = fuse_core();
  std::cout << "Fuzed to core " << core_id << ".\n";

  add_default_counters();

  start_counters();

  // read_begin_counters();
}
// http://researcher.watson.ibm.com/researcher/files/us-ajvega/FastPath_Weaver_Talk.pdf
void PerfStats::add_default_counters() {
  auto PERF_COUNT_RAW_FLOPS = pme_t(0x10, 0x01);
  // auto UOPS_ISSUED_ANY = pme_t(0x0e,0x01);
  // auto UOPS_ISSUED_STALL_CYCLES =pme_t(0x0e,0x01,1,1);

  // auto RESOURCE_STALLS_ANY = pme_t(0xa2,0x01);

  // https://software.intel.com/sites/default/files/managed/a4/60/325384-sdm-vol-3abcd.pdf
  // p 698, vol 3b 18-4
  auto UOPS_RETIRED_ALL = pme_t(0xc0, 0x00);

  auto LLC_REFERENCE = pme_t(0x2e, 0x4f);
  auto LLC_MISSES = pme_t(0x2e, 0x41);
  auto BRANCH_INSTRUCTION_RETIRED = pme_t(0xc4, 0x00);
  auto BRANCH_MISSES_RETIRED = pme_t(0xc5, 0x00);

  auto UNHALTED_REFERENCE_CYCLES = pme_t(0x3c, 0x01);
  auto UNHALTED_CORE_CYCLES = pme_t(0x3c, 0x00);

  // auto FP_ASSIST_ANY= 0x1531ecau;// pme_t(0xca,0x1e);

  /*The time-stamp counter is contained in a 64-bit MSR. The high-order 32 bits
   * of the MSR are loaded into the EDX register, and the low-order 32 bits are
   * loaded into the EAX register. The processor monotonically increments the time-stamp
   * counter MSR every clock cycle and resets it to 0 whenever the processor is reset. */

  add_counter(PERF_TYPE_USER, PERF_TYPE_TICKS, RDTSC, "Ticks");

  // https://software.intel.com/en-us/forums/software-tuning-performance-optimization-platform-monitoring/topic/595214

  /*This event counts the number of instructions that retire execution. For
   * instructions that consist of multiple uops, this event counts the
   * retirement of the last uop of the instruction. The counter continues
   * counting during hardware interrupts, traps, and in-side interrupt handlers.*/
  add_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS, RDPMC, "Instructions", 1 << 30);

  /*The CPU_CLK_UNHALTED.THREAD event counts the number of core
   * cycles while the logical processor is not in a halt state.
   * If there is only one logical processor in a processor core,
   * CPU_CLK_UNHALTED.CORE counts the unhalted cycles of the processor core.
   * The core frequency may change from time to time due to transitions
   * associated with Enhanced IntelSpeedStep Technology or TM2. For this
   * reason this event may have a changing ratio with regards to time.*/
  add_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES, RDPMC, "Actual CPU cycles", (1 << 30) + 1);

  /*This event counts the number of reference cycles at the TSC rate when
   * the core is not in a halt state and not in a TM stop-clock state. The
   * core enters the halt state when it is running the HLT instruction or the
   * MWAIT instruction. This event is not affected by core frequency changes
   * (e.g., P states) but counts at the same frequency as the time stamp counter.
   * This event can approximate elapsed time while the core was not in a halt
   * state and not in a TM stopclock state.*/
  add_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_REF_CPU_CYCLES, RDPMC, "Reference CPU cycles", (1 << 30) + 2);

  // add_counter(PERF_TYPE_RAW, PERF_COUNT_RAW_FLOPS, IOCTL, "Flops");
  // add_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_REFERENCES, IOCTL, "Cache references");
  // add_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES, IOCTL, "Cache misses");
  // add_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_BRANCH_INSTRUCTIONS, IOCTL, "Branch instructions");
  // add_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_BRANCH_MISSES, IOCTL, "Branch misses");
  // add_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_BUS_CYCLES,IOCTL, "BUS cycles");
  // add_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_STALLED_CYCLES_FRONTEND, IOCTL, "Frontend Stalled cycles");
  // add_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_STALLED_CYCLES_BACKEND, IOCTL, "Backend Stalled cycles");

  // add_counter(PERF_TYPE_RAW, FP_ASSIST_ANY, IOCTL, "Flops1");
  // add_counter(PERF_TYPE_RAW,UOPS_ISSUED_ANY, IOCTL, "UOPS issued any");
  // add_counter(PERF_TYPE_RAW,RESOURCE_STALLS_ANY, IOCTL, "Resource stalls any");
  // add_counter(PERF_TYPE_RAW,UOPS_RETIRED_ALL, IOCTL, "UOPS retired all");
  // add_counter(PERF_TYPE_RAW,UOPS_ISSUED_ANY, IOCTL, "UOPS issued any");
  // add_counter(PERF_TYPE_RAW, LLC_REFERENCE, IOCTL, "LLC Reference");
  // add_counter(PERF_TYPE_RAW,LLC_MISSES, IOCTL, "LLC misses");
  // add_counter(PERF_TYPE_RAW,BRANCH_INSTRUCTION_RETIRED, IOCTL, "Branch instruction retired");
  // add_counter(PERF_TYPE_RAW,BRANCH_MISSES_RETIRED, IOCTL, "Branch misses retired");

  // add_counter(PERF_TYPE_RAW, UNHALTED_CORE_CYCLES, IOCTL, "Core cycles");
  // add_counter(PERF_TYPE_RAW, UNHALTED_REFERENCE_CYCLES, IOCTL, "Unhalted reference cycles");
}

void PerfStats::start_counters() {
  for (size_t idx = 0; idx < _counter_count; ++idx) start_hwctr(idx);
}

void PerfStats::stop_counters() {
  for (size_t idx = 0; idx < _counter_count; ++idx) stop_hwctr(idx);
}

void PerfStats::read_begin_counters_lib() {
  long long cnt;

  for (size_t idx = 0; idx < _counter_count; ++idx) {
    if (_conters_dec[idx].method == RDTSC) {
      _begin[idx] = rdtsc_read_ticks_begin();
    } else if (_conters_dec[idx].method == RDPMC) {
      _begin[idx] = rdpmc_read_hwctr(_conters_dec[idx].rdpmc_eax);
    } else if (_conters_dec[idx].method == IOCTL) {
      //   ioctl(_conters_dec[idx].file, PERF_EVENT_IOC_DISABLE, 0);

      if (read(_conters_dec[idx].file, &cnt, 8) == 8)
        _begin[idx] = cnt;
      else
        _begin[idx] = 0;

      //  ioctl(_conters_dec[idx].file, PERF_EVENT_IOC_ENABLE, 0);
    }
  }
}

void PerfStats::read_end_counters_lib() {
  long long cnt;

  for (size_t idx = 0; idx < _counter_count; ++idx) {
    if (_conters_dec[idx].method == RDTSC) {
      _end[idx] = rdtsc_read_ticks_begin();
    } else if (_conters_dec[idx].method == RDPMC) {
      _end[idx] = rdpmc_read_hwctr(_conters_dec[idx].rdpmc_eax);
    } else if (_conters_dec[idx].method == IOCTL) {
      //   ioctl(_conters_dec[idx].file, PERF_EVENT_IOC_DISABLE, 0);

      if (read(_conters_dec[idx].file, &cnt, 8) == 8)
        _end[idx] = cnt;
      else
        _end[idx] = 0;

      //  ioctl(_conters_dec[idx].file, PERF_EVENT_IOC_ENABLE, 0);
    }
    _delta[idx] = _end[idx] - _begin[idx];
  }

  _samples[_sample_count++] = _delta;

  if (_sample_count == MAX_SAMPLES) record_counters();
}
