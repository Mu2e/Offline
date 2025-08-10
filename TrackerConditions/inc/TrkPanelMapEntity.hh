// clang-format off
#ifndef TrackerConditions_TrkPanelMapEntity_hh
#define TrackerConditions_TrkPanelMapEntity_hh
//
// defines cross-mapping of different indices corresponding to a tracker panel
//

// Mu2e includes
#include "Offline/DbTables/inc/TrkPanelMap.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "fhiclcpp/type_traits.h"

// C++ includes
#include <map>
#include <memory>

namespace mu2e {

  class TrkPanelMapEntity : public ProditionsEntity {
  public:
    enum {
      kMaxPlanes = 100,
      kMaxPanels = 600
    };
      
    typedef std::shared_ptr<TrkPanelMapEntity> ptr_t;
    typedef std::shared_ptr<const TrkPanelMapEntity> cptr_t;
    constexpr static const char* cxname = {"TrkPanelMapEntity"};

    TrkPanelMapEntity();
    virtual ~TrkPanelMapEntity() {}

    const TrkPanelMap::Row* panel_map_by_mnid       (int MnID) const {
      return _tpm_by_mnid[MnID];
    };

    const TrkPanelMap::Row* panel_map_by_offline_ind(int Plane, int Panel) const {
      return _tpm_by_offline[Plane][Panel];
    }
    
    const TrkPanelMap::Row* panel_map_by_online_ind (int DtcID, int Link ) const {
      return _tpm_by_online[DtcID][Link];
    }

    void print(std::ostream&) const override;

    void add  (const TrkPanelMap::Row& Tpm);

  private:
                                        // panel map for a given run
    std::vector<TrkPanelMap::Row> _map;
//-----------------------------------------------------------------------------
// assume less than 600 panels, 100 planes
//-----------------------------------------------------------------------------
    const TrkPanelMap::Row*  _tpm_by_mnid   [kMaxPanels];
    const TrkPanelMap::Row*  _tpm_by_offline[kMaxPlanes][6];
    const TrkPanelMap::Row*  _tpm_by_online [kMaxPlanes][6];
  };

}

#endif
