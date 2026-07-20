#ifndef ARITHMETICOPERATIONS_H_
#define ARITHMETICOPERATIONS_H_

#include <stdexcept> // for std::overflow_error
#include <type_traits>

namespace petri
{

// Helper function to add two numbers with overflow checking
template<typename T>
  T addExact (T x, T y)
  {
    if constexpr (std::is_integral_v<T>) {
      T result;
      if (__builtin_add_overflow (x, y, &result)) {
        throw std::overflow_error ("Overflow in addition");
      }
      return result;
    } else {
      return x + y;
    }
  }

// Main function that promotes types and calls the helper
template<typename T1, typename T2>
  auto addExact (T1 x, T2 y) -> typename std::common_type<T1, T2>::type
  {
    using CommonType = typename std::common_type<T1, T2>::type;
    return addExact (static_cast<CommonType> (x), static_cast<CommonType> (y));
  }

// Helper function to multiply two numbers with overflow checking
template<typename T>
  T multiplyExact (T x, T y)
  {
    if constexpr (std::is_integral_v<T>) {
      T result;
      if (__builtin_mul_overflow (x, y, &result)) {
        throw std::overflow_error ("Overflow in multiplication");
      }
      return result;
    } else {
      return x * y;
    }
  }

// Main function that promotes types and calls the helper
template<typename T1, typename T2>
  auto multiplyExact (T1 x, T2 y) -> typename std::common_type<T1, T2>::type
  {
    using CommonType = typename std::common_type<T1, T2>::type;
    return multiplyExact (static_cast<CommonType> (x),
                          static_cast<CommonType> (y));
  }

} // namespace petri

#include <iomanip>

namespace std {
std::ostream& operator<<(std::ostream& out, __uint128_t n) {
  // decimal conversion by repeated division; a 128-bit value has at most 39 digits
  char buf[40];
  size_t pos = sizeof(buf);
  do {
    buf[--pos] = '0' + static_cast<char>(n % 10);
    n /= 10;
  } while (n > 0);
  out.write(buf + pos, sizeof(buf) - pos);
  return out;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> v) {
  out << "[" ;
  bool first = true;
  for (const auto & e : v) {
    if (first) first = false;
    else out << ", ";
    out << e;
  }
  out << "]\n";
  return out;
}

}
#endif /* ARITHMETICOPERATIONS_H_ */
