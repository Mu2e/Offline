#ifndef _perflib_h_
#define _perflib_h_

#include <linux/perf_event.h>
#include "PerfLib/inc/common.hh"
#define MAX_COUNTERS 4
#define MAX_SAMPLES 250
namespace perf {

struct pme_t {
  pme_t(unsigned int event_, unsigned int umask_) : event(event_), umask(umask_) {}
  pme_t(unsigned int event_, unsigned int umask_, unsigned int cmask_) : event(event_), umask(umask_), cmask(cmask_) {}
  pme_t(unsigned int event_, unsigned int umask_, unsigned int cmask_, unsigned int inv_)
      : event(event_), umask(umask_), inv(inv_), cmask(cmask_) {}

  operator unsigned int() const { return _raw; }

  union {
    struct {
      /*This  field  is  used  to  select  the  logic  unit to  detect  the
      *performance  monitoring  event  to  be
      *monitored. So the values to be filled in this field is determined by the
      *architecture.*/
      unsigned int event : 8;

      /*The  logic  unit  selected  by  the  Event  Select  field might  be
      *capable  of  monitoring multiple events.
      *So this UMASK field is used to select one of those events which can be
      *monitored
      *by the logic unit. So based on the logic unit selected the UMASK field
      *may have one fixed
      *value or multiple values, which is dependent on the architecture.*/
      unsigned int umask : 8;

      /*This  flag
      *,  if  set,  tells  the  logic  unit  to  monitor  events  which  happen
      *when  the  processor  is
      *running in the User privilege level i.e. levels 1 through 3.*/
      unsigned int usr : 1;

      /*his  flag,  if  set,  tells  the  logic  unit  to  monitor  events which
      *happen  when  the  processor  is
      *running  in  the  highest  privilege  level
      *i.e.  level  0.  This  flag  and  the  USR  flag  can  be  used
      *together to monitor or count all the events.*/
      unsigned int os : 1;

      /*This flag when set counts the number of times the selected event has
       * gone from low to high state. */
      unsigned int e : 1;

      /*This  flag  when  set  increments  the  counter  and  toggles  the  PMi
       *pin  when  the  monitored event occurs.
       *And when not set, it toggles the PMi pin only when the counter
       *overflows.*/
      unsigned int pc : 1;

      /*When  this  flag  is  set  the  processor  raises  an  interrupt  when
       * the  performance  monitoring counter
       * overflows.*/
      unsigned int intr : 1;

      /*reserved*/
      unsigned int _reserved1 : 1;

      /*This  flag  when  set  enables  the  performance  monitoring  counters
       * for  the  event  and  when clear
       * disables
       * the counters.*/
      unsigned int en : 1;

      /*This  flag  when  set  inverts  the  output  of  CMASK  comparison.
       * This  enables  the  user  to  set both
       * greater than and less than comparisons between CMASK and the counter
       * value.*/
      unsigned int inv : 1;

      /*If this field has a value more tlshan zero then that value is compared
       * to the number of events generated
       * in  one  clock cycle.  If  the  events generated  is  more than the
       * CMASK  value  then the counter is
       * incremented by one else the counter is not incremented.*/
      unsigned int cmask : 8;
    };
    unsigned int _raw : 32;
  };
} __attribute__((packed));

static_assert(sizeof(pme_t) == 4, "Wrong size of pme_t");

struct counter_desc_t {
  unsigned int event_mask;
  unsigned int type;
  int file;
  int method;
  perf_event_attr attr;
  unsigned int rdpmc_eax;
  std::string name;
};

using counters_desc_t = std::array<counter_desc_t, MAX_COUNTERS>;
using counters_t = std::array<unsigned long long, MAX_COUNTERS>;
using samples_t = std::array<counters_t, MAX_SAMPLES>;

using metric_t = std::pair<std::string, double>;
using metrics_t = std::vector<metric_t>;

class PerfStats {
 public:
  PerfStats(std::string const& observation_name) {
    init();
    _observation_name = observation_name;
  }
  ~PerfStats() {
    stop_counters();
    record_counters();
  }

  std::string report_counters();
  void record_counters();

  // void record_counters(std::string const&, size_t);
  void record_metrics(metrics_t const&, std::string const&, unsigned int);

  void read_begin_counters_inlined();
  void read_end_counters_inlined();

  void read_begin_counters_lib();
  void read_end_counters_lib();

 private:
  int fuse_core(int core_id = -1);
  void init();

  void start_hwctr(size_t);
  void stop_hwctr(size_t);
  unsigned long rdpmc_read_hwctr(unsigned int);
  unsigned long long rdtsc_read_ticks_begin();
  unsigned long long rdtsc_read_ticks_end();
  void add_counter(unsigned int, unsigned int, int, std::string const&, unsigned int = 0);
  void add_default_counters();

  void start_counters();
  void stop_counters();

 public:
  size_t _counter_count = 0;
  size_t _sample_count = 0;
  counters_t _begin;
  counters_t _end;
  counters_t _delta;
  samples_t _samples;
  size_t _written_sample_count = 0;
  counters_desc_t _conters_dec;
  std::string _observation_name;
};

inline unsigned long PerfStats::rdpmc_read_hwctr(unsigned int ecx) {
  unsigned int hi, lo;
  __asm__ volatile("rdpmc" : "=a"(lo), "=d"(hi) : "c"(ecx));
  return ((unsigned long)lo) | (((unsigned long)hi) << 32);
}

inline unsigned long long PerfStats::rdtsc_read_ticks_begin() {
  unsigned hi, lo;
  // http://www.intel.com/content/www/us/en/embedded/training/ia-32-ia-64-benchmark-code-execution-paper.html
  __asm__ __volatile__(
      "CPUID\n\t"
      "RDTSC\n\t"
      "mov %%edx, %0\n\t"
      "mov %%eax, %1\n\t"
      : "=r"(hi), "=r"(lo)::"%rax", "%rbx", "%rcx", "%rdx");

  return ((unsigned long long)lo) | (((unsigned long long)hi) << 32);
}

inline unsigned long long PerfStats::rdtsc_read_ticks_end() {
  unsigned hi, lo;
  // http://www.intel.com/content/www/us/en/embedded/training/ia-32-ia-64-benchmark-code-execution-paper.html
  __asm__ __volatile__(
      "RDTSCP\n\t"
      "mov %%edx, %0\n\t"
      "mov %%eax, %1\n\t"
      "CPUID\n\t"
      : "=r"(hi), "=r"(lo)::"%rax", "%rbx", "%rcx", "%rdx");

  return ((unsigned long long)lo) | (((unsigned long long)hi) << 32);
}

enum other_types { PERF_TYPE_USER = PERF_TYPE_MAX + 1, PERF_TYPE_TICKS };
enum methods { IOCTL = 0, RDPMC, RDTSC };

inline void PerfStats::read_begin_counters_inlined() {
  _begin[0] = rdtsc_read_ticks_begin();
  _begin[1] = rdpmc_read_hwctr(_conters_dec[1].rdpmc_eax);
  _begin[2] = rdpmc_read_hwctr(_conters_dec[2].rdpmc_eax);
  _begin[3] = rdpmc_read_hwctr(_conters_dec[3].rdpmc_eax);
}

inline void PerfStats::read_end_counters_inlined() {
  _delta[0] = rdtsc_read_ticks_begin() - _begin[0];
  _delta[1] = rdpmc_read_hwctr(_conters_dec[1].rdpmc_eax) - _begin[1];
  _delta[2] = rdpmc_read_hwctr(_conters_dec[2].rdpmc_eax) - _begin[2];
  _delta[3] = rdpmc_read_hwctr(_conters_dec[3].rdpmc_eax) - _begin[3];

  _samples[_sample_count++] = _delta;

  if (_sample_count == MAX_SAMPLES) record_counters();
}

/*
inline void PerfStats::read_begin_counters() {

  for (size_t idx = 0; idx < _counter_count; ++idx) {
    if (_conters_dec[idx].method == RDTSC) {
      _begin[idx] = rdtsc_read_ticks_begin();
    }else if (_conters_dec[idx].method == RDPMC) {
      _begin[idx] = rdpmc_read_hwctr(_conters_dec[idx].rdpmc_eax);
    }else if (_conters_dec[idx].method == IOCTL) {
      //   ioctl(_conters_dec[idx].file, PERF_EVENT_IOC_DISABLE, 0);
   long long cnt;

      if (read(_conters_dec[idx].file, &cnt, 8) == 8)
        _begin[idx] = cnt;
      else
        _begin[idx] = 0;

      //  ioctl(_conters_dec[idx].file, PERF_EVENT_IOC_ENABLE, 0);
    }
  }
}

inline void PerfStats::read_end_counters() {

  for (size_t idx = 0; idx < _counter_count; ++idx) {
    if (_conters_dec[idx].method == RDTSC) {
      _end[idx] = rdtsc_read_ticks_begin();
    }else if (_conters_dec[idx].method == RDPMC) {
      _end[idx] = rdpmc_read_hwctr(_conters_dec[idx].rdpmc_eax);
    }else if (_conters_dec[idx].method == IOCTL) {
      //   ioctl(_conters_dec[idx].file, PERF_EVENT_IOC_DISABLE, 0);
   long long cnt;

      if (read(_conters_dec[idx].file, &cnt, 8) == 8)
        _end[idx] = cnt;
      else
        _end[idx] = 0;

      //  ioctl(_conters_dec[idx].file, PERF_EVENT_IOC_ENABLE, 0);
    }
     _delta[idx] = _end[idx] - _begin[idx];
  }

  _samples[_sample_count++]=_delta;

   if(_sample_count==MAX_SAMPLES)
     record_counters();

}
*/

}  // namespace perf
#endif /* _perflib_h_ */
