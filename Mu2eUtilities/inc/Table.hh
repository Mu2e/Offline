#ifndef Mu2eUtilities_Table_hh
#define Mu2eUtilities_Table_hh
//
// Free function for loading table two-column table
//
//
// Original author: Kyle Knoepfel

// C++ includes
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

// Framework includes
#include "cetlib_except/exception.h"

namespace mu2e {
  
  // global templates for ease of use
  template<const unsigned M> using Value    = std::array<double,M>;
  template<const unsigned N> using TableRow = std::pair<double,Value<N-1>>;
  template<const unsigned N> using TableVec = std::vector<TableRow<N>>;

  // definition of class
  template <const unsigned N> class Table {
    
    typedef typename TableVec<N>::const_iterator table_const_iterator;

  public:
    
    // Constructors
    Table(){}
    explicit Table( const TableVec<N>& tablevec ) : rawTable_( tablevec ) {}

    // Accessors
    unsigned            getNrows()         const;
    unsigned            getNcols()         const;
    const TableRow<N>&  getRow( const unsigned key ) const;

    const TableVec<N>& rawTable()          const;
    TableVec<N> getShape( const double key1, const double key2, const double res, const double binCorr = 1. ) const;
    TableVec<N> getShape( const std::vector<double>& keys, const double binCorr = 1. ) const;
    
    Value<N-1> getValueAtKey        ( const double key, const double binCorr = 1. ) const;
    unsigned   getLowerBoundRowIndex( const double key ) const;
    
    // Helper utilities
    void printRow( unsigned i ) const;
    void printTable()           const;

    // Modifiers
    void renormalizeShape( const double norm           );
    void interpolateShape( const double key1, const double key2, const double res, const double binCorr = 1. );
    void interpolateShape( const std::vector<double>& keys, const double binCorr = 1. );
    void replaceShape    ( const TableVec<N>& newtable = TableVec<N>(), const bool sortA = true );

    // Operators and ordering functions
    double operator()(unsigned i, unsigned j=1)  const { 
      double val(-1.);
      if ( j == 0 ) val = rawTable_.at(i).first;
      else          val = rawTable_.at(i).second.at(j-1);
      return val;
    }

    struct SortAscend  { bool operator()(const TableRow<N>& row1, const TableRow<N>& row2) { return row1.first < row2.first; } };
    struct SortDescend { bool operator()(const TableRow<N>& row1, const TableRow<N>& row2) { return row1.first > row2.first; } };

    struct CompareKeys { bool operator()(const TableRow<N>& row, const double key){return row.first < key; } };

  private:
    TableVec<N> rawTable_;

    table_const_iterator getLowerBoundRow ( double key ) const; 
    
    template <const unsigned M, const bool SA> 
    friend Table<M> loadTable( const std::string& tableFile );
    
  };
    


  // end of interface
  //====================================================================================
  // generic implementation below
  // explicit specializations in src/Table.cc file

  template <const unsigned N> unsigned Table<N>::getNrows()               const { return rawTable_.size(); }
  template <const unsigned N> unsigned Table<N>::getNcols()               const { return N; }
  template <const unsigned N> const TableRow<N>& Table<N>::getRow(const unsigned i) const { return rawTable_.at(i);}

  template <const unsigned N> const TableVec<N>& Table<N>::rawTable()     const { return rawTable_; }
  

  template <const unsigned N> typename Table<N>::table_const_iterator Table<N>::getLowerBoundRow( double key ) const {
    return std::lower_bound( rawTable_.begin(), rawTable_.end(), key, CompareKeys() );
  }
  
  template <const unsigned N> unsigned Table<N>::getLowerBoundRowIndex( const double key ) const {
    return std::distance( rawTable_.begin(), getLowerBoundRow( key ) );
  }

  template <const unsigned N> void Table<N>::printRow( unsigned i ) const {
    std::cout << rawTable_.at(i).first << " " ;
    for ( auto const & val : rawTable_.at(i).second )
      std::cout << val << " " ;
    std::cout << std::endl;
  }
  
  template <const unsigned N> void Table<N>::printTable() const {
    for ( auto const & it : rawTable_ ) {
      std::cout << it.first << " : " ;
      std::for_each ( it.second.begin(), it.second.end(), [](double val){ std::cout << val << " " ; } );
      std::cout << std::endl;
    }
  }

  template <const unsigned N> void Table<N>::replaceShape( const TableVec<N>& newtable, const bool sortA ) {

    rawTable_.clear();
    if ( newtable.empty() ) return;

    std::copy( newtable.begin(), newtable.end(), rawTable_.begin() );
    if ( sortA ) std::sort( rawTable_.begin(), rawTable_.end(), typename Table<N>::SortAscend()  );
    else         std::sort( rawTable_.begin(), rawTable_.end(), typename Table<N>::SortDescend() );

  }

  //-------------- free function, friend to Table class --------------------------------------------------------
  template <const unsigned N, const bool SA = true>
  Table<N> loadTable( const std::string& tableFile ) {
      
    std::fstream intable(tableFile.c_str(),std::ios::in);
    if ( !intable.is_open() ) {
      throw cet::exception("ProductNotFound")
        << "No Tabulated spectrum table file found";
    }
    
    Table<N> tmp_table;

    // Load table
    while ( !intable.eof() ) {
      TableRow<N> tableRow;
      intable >> tableRow.first;              // Get key first
      std::for_each( tableRow.second.begin(), // Now fill the values
                     tableRow.second.end(), 
                     [&](double& d){
                       intable >> d;
                     } );
      if ( !intable.eof() ) {
        tmp_table.rawTable_.emplace_back( tableRow.first, Value<N-1>( tableRow.second ) );
      }
      
    }
    
    // Optional sorting, with first row corresponding to largest key
    // (i.e. first-column) 
    if ( SA ) std::sort( tmp_table.rawTable_.begin(), tmp_table.rawTable_.end(), typename Table<N>::SortAscend()  );
    else      std::sort( tmp_table.rawTable_.begin(), tmp_table.rawTable_.end(), typename Table<N>::SortDescend() );
    
    return tmp_table;
    
  }
 
} // end of namespace mu2e

#endif /* Mu2eUtilities_Table_hh */
