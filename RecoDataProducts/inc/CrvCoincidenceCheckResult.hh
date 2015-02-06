#ifndef RecoDataProducts_CrvCoincidenceCheckResult_hh
#define RecoDataProducts_CrvCoincidenceCheckResult_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//


namespace mu2e 
{
  class CrvCoincidenceCheckResult
  {
    public:

    CrvCoincidenceCheckResult() {}

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

#endif /* RecoDataProducts_CrvCoincidenceCheckResult_hh */
