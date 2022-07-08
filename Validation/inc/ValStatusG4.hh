
#ifndef ValStatusG4_HH_
#define ValStatusG4_HH_

#include "Offline/MCDataProducts/inc/StatusG4.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValStatusG4 {
 public:
  ValStatusG4(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const StatusG4& obj, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;

  TH1D* _hVer;
  TH1D* _hstat;
  TH1D* _hnTrk;
  TH1D* _hnTrk2;
  TH1D* _hnTrk3;
  TH1D* _hover;
  TH1D* _hkill;
  TH1D* _hkillfp;
  TH1D* _hCPU1;
  TH1D* _hCPU2;
  TH1D* _hCPU3;
  TH1D* _hWall1;
  TH1D* _hWall2;
  TH1D* _hWall3;
};
}  // namespace mu2e

#endif
