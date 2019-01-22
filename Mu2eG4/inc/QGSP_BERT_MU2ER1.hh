// Physics list inheriting from QGSP_BERT using regions etc.
// Uses QGSP_BERT
// Original author K.L. Genser October 2018
//
#ifndef Mu2eG4_QGSP_BERT_MU2ER1_h
#define Mu2eG4_QGSP_BERT_MU2ER1_h 1

#include "QGSP_BERT.hh"

namespace mu2e {

  class QGSP_BERT_MU2ER1: public QGSP_BERT // could try to templetize it here
  {
  public:

    explicit QGSP_BERT_MU2ER1( const fhicl::ParameterSet& pSet, G4int verbose = 1 );

    virtual ~QGSP_BERT_MU2ER1(){};

    virtual void SetCuts() override;

  private:

    const fhicl::ParameterSet& pset;

  };

}
#include "QGSP_BERT_MU2ER1.icc"

#endif
