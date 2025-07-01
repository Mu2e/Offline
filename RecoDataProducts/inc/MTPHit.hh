#ifndef RecoDataProducts_MTPHit_hh
#define RecoDataProducts_MTPHit_hh

#include <vector>

namespace mu2e
{

  class MTPHit
  {

  public:

    // constructors
    MTPHit() : _time(0), _channelID(0) {};
    MTPHit(float Time, int ChannelID) : _time(Time), _channelID(ChannelID) {};

    // accessors
    float const& time()	      const { return _time;	 }
    int	  const& channelID()  const { return _channelID; }

  private:

    // data members
    float _time;      // time of hit
    int	  _channelID; // paddle

  };

  typedef std::vector<mu2e::MTPHit> MTPHitCollection;

} // namespace mu2e

#endif /* RecoDataProducts_MTPHit_hh */
