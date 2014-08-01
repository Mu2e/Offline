#ifndef CalorimeterGeom_DiskCrystalPosUtil_hh
#define CalorimeterGeom_DiskCrystalPosUtil_hh
//
// $Id: DiskCrystalPosUtil.hh,v 1.2 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//

//C++ includes
#include <vector>
#include <map>


namespace mu2e {


    class DiskCrystalPosUtil {


	public:

            DiskCrystalPosUtil(){}

	    void Fill(int iCrystal); 
	    void Fill(int iCrystal, int iMap, int l, int k); 
	    
	    int mapToCrystal(int iMap)     const;    
	    int crystalToMap(int iCrystal) const;
	    int crystalToL(int iCrystal)   const; 
	    int crystalToK(int iCrystal)   const; 
	    int LKToCrystal(int l, int k)  const;
	                                  

	private:

            std::vector<int> _mapToCrystal;
	    std::vector<int> _crystalToMap;
	    std::vector<int> _crystalToL;
	    std::vector<int> _crystalToK;
	    std::map< int, int> _LKToCrystal;	    
    };

}
#endif

