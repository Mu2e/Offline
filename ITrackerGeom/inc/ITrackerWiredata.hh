#ifndef ITrackerGeom_ITrackerWiredata_hh
#define ITrackerGeom_ITrackerWiredata_hh

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  DCH wireposition class                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObjArray.h"

class ITrackerWiredata {

  public:

    ITrackerWiredata();

    virtual ~ITrackerWiredata();

    TObjArray *PosMatrix; //PosMatrix
    Int_t     NcelLayer;  //NcelLayer
    Float_t   *epsilon;   //[NcelLayer]  epsilon
    Float_t   *alfa;      //[NcelLayer]  alfa
    Float_t   *radius_z0; //[NcelLayer]  layer radius at z=0

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

#endif /* ITrackerGeom_ITrackerWiredata_hh */
