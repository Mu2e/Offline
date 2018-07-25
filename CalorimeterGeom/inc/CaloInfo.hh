//
// Contains data for the calorimeter class. The non-critical fields are saved into maps, and a few performance critical fields
// accessed throughout the code are cached for efficiency.
//
// Original author B. Echenard
//
#ifndef CalorimeterGeom_CaloInfo_hh
#define CalorimeterGeom_CaloInfo_hh

#include "cetlib_except/exception.h"

#include <vector>
#include <map>
#include <string>


namespace mu2e {
   
    //helper class
    template <typename T> class CaloInfoData
    {
       public:     
          CaloInfoData() : data_() {};
          ~CaloInfoData() {};

          const T& get(const std::string& key) const
          { 
             auto iter = data_.find(key); 
             if (iter == data_.end()) throw cet::exception("CaloInfo") << " unknown element "<<key<<"\n";
             return iter->second;
          };

          void set(const std::string key, const T& value) {data_[key] = value;} 


       private:    
          std::map<std::string,T> data_;
    };

    
    
    
    
    
    
    class CaloInfo {

       public:

           CaloInfo() : dataBool_(),dataInt_(),dataDouble_(),dataVDouble_(),nROPerCrystal_(0)
           {}
	     
           ~CaloInfo() {}
           
	   void nROPerCrystal(int value)             {nROPerCrystal_ = value;}
  	   int  nROPerCrystal()                const {return nROPerCrystal_;}         
	   int  crystalByRO(int roid)          const {return (roid/nROPerCrystal_);}
	   int  ROBaseByCrystal(int crystalId) const {return (crystalId*nROPerCrystal_);}

           void set(const std::string& key, int value)                       {dataInt_.set(key,value);}
           void set(const std::string& key, double value)                    {dataDouble_.set(key,value);}
           void set(const std::string& key, const std::vector<double>& value){dataVDouble_.set(key,value);}
           
           const bool                getBool(const std::string& key)    const {return dataBool_.get(key);}
           const int                 getInt(const std::string& key)     const {return dataInt_.get(key);}
           const double              getDouble(const std::string& key)  const {return dataDouble_.get(key);}
           const std::vector<double> getVDouble(const std::string& key) const {return dataVDouble_.get(key);}
          
           //helper function
           double crystalVolume() const {return dataDouble_.get("crystalXYLength")*dataDouble_.get("crystalXYLength")*dataDouble_.get("crystalZLength");}


       private:

          CaloInfoData<bool>                dataBool_;
          CaloInfoData<int>                 dataInt_;
          CaloInfoData<double>              dataDouble_;
          CaloInfoData<std::vector<double>> dataVDouble_;
          int                               nROPerCrystal_;
          
          
          
                   
     };

}    

#endif 
