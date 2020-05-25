#ifndef CaloPulseCache_HH
#define CaloPulseCache_HH

#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include <vector>
#include <string> 


namespace mu2e {


   class CaloPulseCache {

       public:

          CaloPulseCache();
          ~CaloPulseCache() {};

	  void   initialize();
          double evaluate(double x);

          const std::vector<double>&   cache()      const {return cache_;}
          double                       cache(int i) const {return cache_.at(i);}
          double                       cacheSize()  const {return cacheSize_;}
          double                       deltaT()     const {return deltaT_;}
          double                       step()       const {return step_;}


      private:

	 std::vector<double> cache_;
	 double              cacheSize_;
	 double              deltaT_;
	 double              step_;
   };

}
#endif
