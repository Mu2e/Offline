#ifndef RecoDataProducts_CRVCoincidenceCheckResult_hh
#define RecoDataProducts_CRVCoincidenceCheckResult_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//


namespace mu2e 
{
  class CRVCoincidenceCheckResult
  {
    public:

    CRVCoincidenceCheckResult() {}

    void SetCoincidence(bool coincidence)
    {
      _coincidence = coincidence;
    }

    const bool GetCoincidence() const 
    {
      return _coincidence;
    }

    private:

    bool   _coincidence;
  };
}

#endif /* RecoDataProducts_CRVCoincidenceCheckResult_hh */
