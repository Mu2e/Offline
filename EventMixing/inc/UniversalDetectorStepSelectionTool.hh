// Ed Callaghan
// Select all detector steps
// November 2024

#ifndef EventMixing_UniversalDetectorStepSelectionTool_hh
#define EventMixing_UniversalDetectorStepSelectionTool_hh

// art
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/EventMixing/inc/DetectorStepSelectionTool.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"

namespace mu2e{
  class UniversalDetectorStepSelectionTool: public DetectorStepSelectionTool{
    public:
      struct Config{
        /**/
      };

      using Parameters = art::ToolConfigTable<Config>;
      UniversalDetectorStepSelectionTool(const Parameters&);
      UniversalDetectorStepSelectionTool() = default;
     ~UniversalDetectorStepSelectionTool() = default;

      virtual bool Select(const CaloShowerStep&) override final;
      virtual bool Select(const CrvStep&)        override final;
      virtual bool Select(const StrawGasStep&)   override final;

    protected:
      template<typename T>
      bool select(const T&);

    public:
      /**/
  };

  // trivial filter; select everything
  template<typename T>
  bool UniversalDetectorStepSelectionTool::select(const T& step){
    bool rv = true;
    return rv;
  }

} // namespace mu2e

#endif
