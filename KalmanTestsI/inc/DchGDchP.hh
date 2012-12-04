#ifndef DCHGDCHP_HH
#define DCHGDCHP_HH

#include "DchGeomBase/DchGDch.hh"
#include <string>

class DchGDchP: public DchGDch{
public:
  DchGDchP(std::string materialfile="", bool doDetailedWrSupp=true, bool useSingleCellElemDescr=false, bool gas_wire_extWght=false);

private:
  void makeDetailedWrSupp(std::vector<DchPhiSegmCyl> &epCylSubs,
                          std::vector<double*> &ePPosSubs,
                          std::vector<std::string> &ePMaterialSubs,
                          double minZ=0.0, double totWidth=0.0, double rMin=0.0, double rMax=0.0, bool isFEP=false);

  void makeInEPCables(std::vector<DchPhiSegmCyl> &epCylSubs,
                          std::vector<double*> &ePPosSubs,
                          std::vector<std::string> &ePMaterialSubs,
                          double minZ=0.0, double totWidth=0.0, double rMin=0.0, double rMax=0.0, bool isREP=false);

};

#endif
