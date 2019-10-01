// simple counting histogram class based on std::Array.
// Count should be a unsigned integer type (uint16_t, ...)
//  Original author: Dave Brown (LBNL) Feb 2019
#ifndef DataProducts_AHist_HH
#define DataProducts_AHist_HH
namespace mu2e {
  template <class Count, std::size_t nbins> class AHist{
    public:
      AHist(float low, float hi) : _hlow(low), _hbin((hi-low)/nbins), _hist{0} {}
      float binSize() const { return _hbin;}
      float binLowEdge(size_t ibin) const { return _hlow + binSize()*ibin; }
      float binMid(size_t ibin) const { return _hlow + binSize()*(float(ibin) + 0.5); }
      float binHighEdge(size_t ibin) const { return binLowEdge(ibin+1); }
      int binIndex(float value) const { return (int)floorf( (value-_hlow)/binSize()); }
      Count binContents(size_t ibin) const { return ibin < nbins ? _hist[ibin] : 0; }
      Count binContents(float value) const {
	auto ibin = binIndex(value);
	if(ibin >= 0 && ibin < nbins) return _hist[Count(ibin)];
      }
      void fill(float value,Count increment=1) { 
	auto ibin = binIndex(value);
	if(ibin >= 0 && ibin < (int)nbins)_hist[Count(ibin)]+= increment;
      }
    private:
      float _hlow, _hbin; // low edge and bin size
      std::array<Count,nbins> _hist; // histogram counts
  }; 
}
#endif 
