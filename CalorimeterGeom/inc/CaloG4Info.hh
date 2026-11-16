//
// Contains data for the calorimeter class. The non-critical fields are stored in a single keyed
// store; a few performance-critical fields accessed throughout the code are cached separately.
//
// Access is typed: set<T>/get<T> with T one of the supported value types (see Value below).
// See wiki for list of valid keys.
//
// Original author B. Echenard
//
#ifndef CalorimeterGeom_CaloG4Info_hh
#define CalorimeterGeom_CaloG4Info_hh

#include "cetlib_except/exception.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>


namespace mu2e {

  class CaloG4Info {
    public:
      CaloG4Info() = default;

      template <typename T>
      void set(const std::string& key, T value)
      {
          static_assert(isSupported<T>, "CaloG4Info::set - unsupported value type");
          data_[key] = std::move(value);
      }

      template <typename T>
      const T& get(const std::string& key) const
      {
          static_assert(isSupported<T>, "CaloG4Info::get - unsupported value type");
          const auto iter = data_.find(key);
          if (iter == data_.end())
             throw cet::exception("CaloG4Info") << " unknown element " << key << "\n";
          if (!std::holds_alternative<T>(iter->second))
             throw cet::exception("CaloG4Info") << " element " << key << " requested with the wrong type\n";
          return std::get<T>(iter->second);
      }


    private:
      using Value = std::variant<bool, int, double,
                                 std::vector<int>, std::vector<double>,
                                 std::string, CLHEP::Hep3Vector>;

      template <typename T, typename V> struct IsAlt : std::false_type {};
      template <typename T, typename... Ts>
      struct IsAlt<T, std::variant<Ts...>> : std::disjunction<std::is_same<T,Ts>...> {};

      template <typename T> static constexpr bool isSupported = IsAlt<T,Value>::value;

      std::map<std::string,Value> data_;
  };
}

#endif
