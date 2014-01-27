#ifndef Mu2eUtilities_Table_hh
#define Mu2eUtilities_Table_hh
//
// Free function for loading table two-column table
//
// $Id: Table.hh,v 1.3 2014/01/27 22:20:17 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/01/27 22:20:17 $
//
// Original author: Kyle Knoepfel

// Mu2e includes
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

// C++ includes
#include <array>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

// Framework includes
#include "cetlib/exception.h"

namespace mu2e {
  
  // global template alias for readability
  template <const unsigned N>
  using TableRow = std::array<double,N>;
  
  // definition of class
  template <const unsigned N> 
  class Table {
    
    typedef typename std::vector<TableRow<N>> RawTable;
    typedef typename RawTable::const_iterator table_const_iterator;

  public:

    unsigned            getNrows()        const { return _rawTable.size(); }
    unsigned            getNcols()        const { return N; }
    const TableRow<N>&  row( unsigned i ) const { return _rawTable.at(i);  }
    const RawTable&     rawTable()        const { return _rawTable;     }

    int findLowerBoundRow( double e ) const {
      auto const & it = getLowerBoundRow( e );
      return ( it - _rawTable.begin() ); 
    }

    double operator()(unsigned i, unsigned j=1)         const { return _rawTable.at(i).at(j); }

    void   printRow( unsigned i )                       const {
      for ( auto const & val : _rawTable[i] )
        std::cout << val << " " ;
      std::cout << std::endl;
    }

    void   printTable()                               const {
      for ( auto const & it : _rawTable ) {
        std::for_each ( it.begin(), it.end(), [](double val){ std::cout << val << " " ; } );
        std::cout << std::endl;
      }
    }

  private:
    RawTable _rawTable;

    const table_const_iterator getLowerBoundRow( double e ) const {
      return std::lower_bound( _rawTable.begin(), _rawTable.end(), e,
                               [](TableRow<2> it, double E){ return it[0] > E; } );
    }
    
    template <const unsigned M> 
    friend Table<M> loadTable( const std::string& tableFile, const bool sort );
    
  };
    
  // Free function to load table
  template <const unsigned N>
  Table<N> loadTable( const std::string& tableFile, 
                      const bool sort = true ) {
      
    ConfigFileLookupPolicy findConfig;
    std::string filename = findConfig( tableFile.c_str() );
      
    std::fstream intable(filename.c_str(),std::ios::in);
    if ( !intable.is_open() ) {
      throw cet::exception("ProductNotFound")
        << "No Tabulated spectrum table file found";
    }
    
    Table<N> tmp_table;

    // Load table
    while ( !intable.eof() ) {
      TableRow<N> tableRow;
      std::for_each( tableRow.begin(), 
                     tableRow.end(), 
                     [&](double& d){
                       intable >> d;
                     } );
      if ( !intable.eof() )
        tmp_table._rawTable.push_back( tableRow );
    }

    // Optional sorting, with first row corresponding to largest key
    // (i.e. first-column) value energy (to be compatible with
    // existing code)
    if ( sort ) {
      std::sort( tmp_table._rawTable.begin(), 
                 tmp_table._rawTable.end(),
                 [](TableRow<N> lhs, TableRow<N> rhs) {
                   return lhs[0] > rhs[0]; 
                 } );
    }

    return tmp_table;

  }

} // end of namespace mu2e

#endif /* Mu2eUtilities_Table_hh */
