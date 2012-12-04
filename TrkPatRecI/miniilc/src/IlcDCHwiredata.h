#ifndef ILCDCHWIREDATA_H
#define ILCDCHWIREDATA_H

#include <TObject.h>
 
class IlcDCHwiredata : public TObject { 

  public: 

    IlcDCHwiredata(){
      PosMatrix = new TObjArray(0);
      NcelLayer = 0;
      epsilon = 0x0;
      alfa = 0x0;
      radius_z0 = 0x0;
    }

    virtual ~IlcDCHwiredata(){
      PosMatrix->Delete();
      
      NcelLayer = 0;
      if (epsilon) delete [] epsilon;
      if (alfa) delete [] alfa;
      if (radius_z0) delete [] radius_z0;
    }
  
    TObjArray *PosMatrix; //PosMatrix
    int       NcelLayer;  //NcelLayer
    float   *epsilon;   //[NcelLayer]  epsilon
    float   *alfa;      //[NcelLayer]  alfa
    float   *radius_z0; //[NcelLayer]  layer radius at z=0

  ClassDef(IlcDCHwiredata,1)                                // DCH geometry class  
};
 
#endif
