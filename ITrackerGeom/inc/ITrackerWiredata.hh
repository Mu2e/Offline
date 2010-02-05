#ifndef ITRACKERWIREDATA_HH
#define ITRACKERWIREDATA_HH

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  DCH wireposition class                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TObjArray.h"

class ITrackerWiredata : public TObject {

  public:

    ITrackerWiredata();

    virtual ~ITrackerWiredata();

    TObjArray *PosMatrix; //PosMatrix
    Int_t     NcelLayer;  //NcelLayer
    Float_t   *epsilon;   //[NcelLayer]  epsilon
    Float_t   *alfa;      //[NcelLayer]  alfa
    Float_t   *radius_z0; //[NcelLayer]  layer radius at z=0

//  ClassDef(ITrackerWiredata,1)                                // DCH geometry class
};

#endif

#ifndef ILCDCHWIREDATA_HH
#define ILCDCHWIREDATA_HH

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  DCH wireposition class                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class IlcDCHwiredata : public ITrackerWiredata {

  public:

    IlcDCHwiredata():ITrackerWiredata(){}

    virtual ~IlcDCHwiredata(){}
};

#endif
