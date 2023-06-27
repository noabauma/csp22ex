/* Computational Statistical Physics 2021
 * Exercise 8
 * Author: Pascal Engeler, engelerp@phys.ethz.ch
 */
#ifndef VEC3_HPP_INCLUDED
#define VEC3_HPP_INCLUDED
#if defined (_MSC_VER)
#define REALLY_INLINE
#else
#define REALLY_INLINE __attribute__((always_inline))
#endif
#include <cmath>

template<typename element_t>
class Vec3{
public:
  Vec3(): x_(0), y_(0), z_(0) {}
  Vec3(element_t x, element_t y, element_t z): x_(x), y_(y), z_(z) {}
  Vec3(const Vec3<element_t>&) = default;
  Vec3<element_t>& operator=(const Vec3<element_t>&) = default;

  Vec3<element_t>& operator+=(const Vec3<element_t>& rhs) noexcept{
    x_ += rhs.x_;
    y_ += rhs.y_;
    z_ += rhs.z_;
    return *this;
  }
  Vec3<element_t>& operator-=(const Vec3<element_t>& rhs) noexcept{
    x_ -= rhs.x_;
    y_ -= rhs.y_;
    z_ -= rhs.z_;
    return *this;
  }
  Vec3<element_t>& operator*=(const element_t& scalar) noexcept{
    x_ *= scalar;
    y_ *= scalar;
    z_ *= scalar;
    return *this;
  }

  inline const element_t rsq() const noexcept REALLY_INLINE{
    return x_*x_ + y_*y_ + z_*z_;
  }
  inline const element_t r() const noexcept REALLY_INLINE{
    return std::sqrt(rsq());
  }

  inline const element_t x() const noexcept REALLY_INLINE {return x_;}
  inline const element_t y() const noexcept REALLY_INLINE {return y_;}
  inline const element_t z() const noexcept REALLY_INLINE {return z_;}
  inline element_t& x() noexcept REALLY_INLINE {return x_;}
  inline element_t& y() noexcept REALLY_INLINE {return y_;}
  inline element_t& z() noexcept REALLY_INLINE {return z_;}
private:
  element_t x_, y_, z_;
};


template<typename element_t>
Vec3<element_t> operator+(Vec3<element_t> lhs, const Vec3<element_t>& rhs) noexcept{
  lhs += rhs;
  return lhs;
}

template<typename element_t>
Vec3<element_t> operator-(Vec3<element_t> lhs, const Vec3<element_t>& rhs) noexcept{
  lhs -= rhs;
  return lhs;
}

template<typename scalar_t, typename element_t>
Vec3<element_t> operator*(const scalar_t scalar, Vec3<element_t> rhs) noexcept{
  rhs *= element_t(scalar);
  return rhs;
}

template<typename scalar_t, typename element_t>
Vec3<element_t> operator*(Vec3<element_t> lhs, const scalar_t scalar) noexcept{
  lhs *= element_t(scalar);
  return lhs;
}
/*
template<typename element_t>
Vec3<element_t> operator*(const element_t scalar, const Vec3<element_t> rhs) noexcept {
  Vec3<element_t> retvec(rhs);
  retvec *= scalar;
  return retvec;
}

template<typename element_t>
Vec3<element_t> operator*(const Vec3<element_t> lhs, const element_t scalar) noexcept {
  Vec3<element_t> retvec(lhs);
  retvec *= scalar;
  return retvec;
}
*/
#endif
