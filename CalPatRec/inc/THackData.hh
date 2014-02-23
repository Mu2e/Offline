//
#ifndef THackData_hh
#define THackData_hh

#include "TNamed.h"
#include "TH1F.h"

class HelixTraj;
//class TH1F;

class THackData: public TNamed {
public:

  double  fClusterT0;
  double  fClusterZ;
  
  const HelixTraj*  fSHelix;

  double  fData[100];

  //  TH1F *fHdist;

  THackData();
  THackData(const char* Name, const char* Title);
  ~THackData();

  double ClusterT0() { return fClusterT0; }
  double ClusterZ()  { return fClusterZ;  }
  double ClusterX()  { return fData[0];   }
  double ClusterY()  { return fData[1];   }
  int    SeedIndex() { return fData[2];   }
  int    Seed2Index(){ return fData[3];   }
 
  void   SetClusterT0(double T0) { fClusterT0 = T0; }
  void   SetClusterZ (double z)  { fClusterZ  = z;  }
};

#endif
