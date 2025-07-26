#ifndef MCDataProducts_ScorerSummary_hh
#define MCDataProducts_ScorerSummary_hh

//
// Information about scorers managed by Geant4
//
// Original author Bertrand Echenard
//

#include <vector>

namespace mu2e {

  struct ScorerSummary {

    // A default c'tor is required for ROOT.
    ScorerSummary():
      ix_(0u),
      iy_(0u),
      iz_(0u),
      entries_(0u),
      total_(0.0),
      totalSqr_(0.0)
    {}

    ScorerSummary(unsigned ix, unsigned iy, unsigned iz, unsigned entries, double total, double totalSqr):
      ix_(ix),
      iy_(iy),
      iz_(iz),
      entries_(entries),
      total_(total),
      totalSqr_(totalSqr)
    {}

    unsigned ix()         const {return ix_;}
    unsigned iy()         const {return iy_;}
    unsigned iz()         const {return iz_;}
    unsigned entries()    const {return entries_;}
    float    total()      const {return total_;}
    float    totalSqr()   const {return totalSqr_;}

  private:
    unsigned ix_;
    unsigned iy_;
    unsigned iz_;
    unsigned entries_;
    float    total_;
    float    totalSqr_;
  };

  typedef std::vector<mu2e::ScorerSummary> ScorerSummaryCollection;
}

#endif
