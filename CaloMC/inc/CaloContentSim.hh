#ifndef CaloCluster_CaloContentSim_HH_
#define CaloCluster_CaloContentSim_HH_

#include <algorithm>

namespace mu2e {


     class CaloContentSim {

	   public:

	       CaloContentSim(double edep, double time, double mom): edep_(edep), time_(time),mom_(mom) {};
	       
	       void update(double edep, double time, double mom) {edep_ += edep; time_ = std::min(time, time_); mom_ = std::max(mom, mom_);}

	       double edep() const {return edep_;}
	       double time() const {return time_;}
	       double mom()  const {return mom_;}


	   private:
	   
	       double edep_;
	       double time_;
	       double mom_;
     };
} 

#endif
