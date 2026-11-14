#ifndef CRVDQM_inc_CRVDigiDQMConfig_hh
#define CRVDQM_inc_CRVDigiDQMConfig_hh
// Validated FHiCL for CRVDigiDQM: the generic `hists` table and the configured
// digi-time pairs. Offline modules nest it with fhicl::Table; online modules
// build the same table from their ParameterSet.
//
//   dqm : {
//     hists : { ... }                                   # DQMHelpers/README.md
//     dtPairs : [
//       { rocA : 1  febA : 3  fpgaA : 0  fpgaB : 1 },   # one FEB, FPGA 1 - FPGA 0
//       { rocA : 1  febA : 3  fpgaA : 2  fpgaB : 2 },   # one FPGA, its different channels
//       { rocA : 1  febA : 3  fpgaA : 0  rocB : 2  febB : 5  fpgaB : 0 }   # across ROCs
//     ]
//   }
//
// rocB / febB default to rocA / febA. FEB is the ROC port, 1-based; FPGA is
// 0-3. CRVDigiDQM refuses two FEBs of one ROC (see CRVDigiDQM.hh).
//
// Original Author: R. Mina

#include "Offline/CRVDQM/inc/CRVDigiDQM.hh"
#include "Offline/DQMHelpers/inc/DQMHistSetConfig.hh"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Table.h"

#include <vector>

namespace mu2e {

struct CRVDtPairFhicl {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;
  fhicl::Atom<int> rocA{Name("rocA"), Comment("ROC of A, 1-based")};
  fhicl::Atom<int> febA{Name("febA"), Comment("FEB (ROC port) of A, 1-based")};
  fhicl::Atom<int> fpgaA{Name("fpgaA"), Comment("FPGA of A, 0-3")};
  fhicl::OptionalAtom<int> rocB{Name("rocB"), Comment("ROC of B; default rocA")};
  fhicl::OptionalAtom<int> febB{Name("febB"), Comment("FEB (ROC port) of B; default febA")};
  fhicl::Atom<int> fpgaB{Name("fpgaB"), Comment("FPGA of B, 0-3")};
};

struct CRVDigiDQMFhicl {
  fhicl::Table<DQMHistSetFhicl> hists{
      fhicl::Name("hists"),
      fhicl::Comment("copies and publishing; see Offline/DQMHelpers/README.md")};
  fhicl::OptionalSequence<fhicl::Table<CRVDtPairFhicl>> dtPairs{
      fhicl::Name("dtPairs"),
      fhicl::Comment("(ROC, FEB, FPGA) pairs to book a dtPair_* histogram for, t_B - t_A over "
                     "every digi pair; omit for none")};
};

inline std::vector<CRVDigiDQM::DtPair> toDtPairs(const CRVDigiDQMFhicl& c)
{
  std::vector<CRVDigiDQM::DtPair> out;
  std::vector<CRVDtPairFhicl> pairs;
  if (!c.dtPairs(pairs)) {
    return out;
  }
  for (const CRVDtPairFhicl& p : pairs) {
    CRVDigiDQM::DtPair d;
    d.a = {p.rocA(), p.febA(), p.fpgaA()};
    d.b = {p.rocB().value_or(p.rocA()), p.febB().value_or(p.febA()), p.fpgaB()};
    out.push_back(d);
  }
  return out;
}
// Online modules validate their ParameterSet the same way:
//   const fhicl::Table<CRVDigiDQMFhicl> t(ps.get<fhicl::ParameterSet>("digi", {}));
//   CRVDigiDQM dqm(toConfig(t().hists()), configuration, toDtPairs(t()));

} // namespace mu2e

#endif /* CRVDQM_inc_CRVDigiDQMConfig_hh */
