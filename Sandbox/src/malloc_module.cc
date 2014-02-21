// A module to allocate a controlled amount of memory.
//
// Andrei Gaponenko, 2014

#include <string>
#include <iostream>
#include <cstdlib>
#include <cstddef>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

namespace mu2e {

  //================================================================
  class malloc : public art::EDAnalyzer {
  public:
    explicit malloc(fhicl::ParameterSet const& pset);
    ~malloc();
    void analyze(const art::Event& evt) override;
  private:
    void *ptr_;
  };

  //================================================================
  malloc::malloc(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , ptr_()
  {
    size_t nbytes = pset.get<size_t>("nbytes");
    std::cerr<<"malloc_module: allocating "<<nbytes<<" bytes"<<std::endl;
    ptr_ = ::malloc(nbytes);
    std::cerr<<"malloc_module: allocated "<<nbytes<<" bytes: ptr = "<<ptr_<<std::endl;
  }

  //================================================================
  malloc::~malloc() {
    try {
      std::cerr<<"malloc_module: releasing ptr="<<ptr_<<std::endl;
      free(ptr_);
    }
    catch(...) {}
  }

  //================================================================
  void malloc::analyze(const art::Event& event) {
    sleep(1);
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::malloc);
