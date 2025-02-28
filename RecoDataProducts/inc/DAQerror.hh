#ifndef RecoDataProducts_DAQerror_hh
#define RecoDataProducts_DAQerror_hh
//
//
// Contact person Ralf Ehrlich
//

#include "Offline/GeneralUtilities/inc/EnumToStringSparse.hh"
#include <vector>

namespace mu2e
{
  class DAQerrorCodeDetail
  {
    public:

    enum enum_type{unknown=0, byteCountMismatch=1};
    static std::string const& typeName();
    static std::map<enum_type,std::string> const& names();
  };
  typedef EnumToStringSparse<DAQerrorCodeDetail> DAQerrorCode;

  class DAQerror
  {
    public:

    DAQerror() :
               _errorCode(), _fragmentIndex(0) {}

    DAQerror(DAQerrorCode::type errorCode, size_t fragmentIndex) :
               _errorCode(errorCode), _fragmentIndex(fragmentIndex) {}

    DAQerrorCode::type GetErrorCode() const     {return _errorCode;}
    size_t             GetFragmentIndex() const {return _fragmentIndex;}

    private:

    DAQerrorCode::type _errorCode;
    size_t             _fragmentIndex;
  };
  typedef std::vector<mu2e::DAQerror> DAQerrorCollection;
}

#endif /* RecoDataProducts_DAQerror_hh */
