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
  double TheoImode (){ return fData[4];   }
  double TheoRadius(){ return fData[5];   }
  double TheoPhi0  (){ return fData[6];   }
  double TheoTanL  (){ return fData[7];   }
  double Theodfdz  (){ return fData[8];   }  
  int    mode0Points(){ return fData[9];   }  
  double shDz      (){ return fData[10];  }  
  int    goodPoints(){ return fData[11];  }  
  double chi2      (){ return fData[12];  }

  
  void   SetClusterT0(double T0) { fClusterT0 = T0; }
  void   SetClusterZ (double z)  { fClusterZ  = z;  }
};

#endif
