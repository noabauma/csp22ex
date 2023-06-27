/* Computational Statistical Physics 2022
 * Exercise 12
 * Author: Pascal Engeler, engelerp@phys.ethz.ch
 */
#ifndef VEC2_HPP_INCLUDED
#define VEC2_HPP_INCLUDED
#if defined (_MSC_VER)
#define REALLY_INLINE
#else
#define REALLY_INLINE __attribute__((always_inline))
#endif
#include <cmath>

template<typename element_t>
class Vec2{
public:
  Vec2(): x_(0), y_(0) {}
  Vec2(element_t x, element_t y): x_(x), y_(y) {}
  Vec2(const Vec2<element_t>&) = default;
  Vec2<element_t>& operator=(const Vec2<element_t>&) = default;

  Vec2<element_t>& operator+=(const Vec2<element_t>& rhs) noexcept{
    x_ += rhs.x_;
    y_ += rhs.y_;
    return *this;
  }
  Vec2<element_t>& operator-=(const Vec2<element_t>& rhs) noexcept{
    x_ -= rhs.x_;
    y_ -= rhs.y_;
    return *this;
  }
  Vec2<element_t>& operator*=(const element_t& scalar) noexcept{
    x_ *= scalar;
    y_ *= scalar;
    return *this;
  }
  Vec2<element_t>& operator/=(const element_t& scalar){
    x_ /= scalar;
    y_ /= scalar;
    return *this;
  }
  inline Vec2<element_t> operator-() const noexcept REALLY_INLINE{
    return {-x_, -y_};
  }
  inline const element_t rsq() const noexcept REALLY_INLINE{
    return x_*x_ + y_*y_;
  }
  inline const element_t r() const noexcept REALLY_INLINE{
    return std::sqrt(rsq());
  }

  inline const element_t x() const noexcept REALLY_INLINE {return x_;}
  inline const element_t y() const noexcept REALLY_INLINE {return y_;}
  inline element_t& x() noexcept REALLY_INLINE {return x_;}
  inline element_t& y() noexcept REALLY_INLINE {return y_;}
private:
  element_t x_, y_;
};


template<typename element_t>
Vec2<element_t> operator+(Vec2<element_t> lhs, const Vec2<element_t>& rhs) noexcept{
  lhs += rhs;
  return lhs;
}

template<typename element_t>
Vec2<element_t> operator-(Vec2<element_t> lhs, const Vec2<element_t>& rhs) noexcept{
  lhs -= rhs;
  return lhs;
}

template<typename scalar_t, typename element_t>
Vec2<element_t> operator*(const scalar_t scalar, Vec2<element_t> rhs) noexcept{
  rhs *= element_t(scalar);
  return rhs;
}

template<typename scalar_t, typename element_t>
Vec2<element_t> operator*(Vec2<element_t> lhs, const scalar_t scalar) noexcept{
  lhs *= element_t(scalar);
  return lhs;
}

template<typename scalar_t, typename element_t>
Vec2<element_t> operator/(Vec2<element_t> lhs, const scalar_t scalar){
  lhs /= element_t(scalar);
  return lhs;
}

template<typename element_t>
element_t dot(const Vec2<element_t>& lhs, const Vec2<element_t>& rhs) noexcept{
  return lhs.x() * rhs.x() + lhs.y() * rhs.y();
}
#endif
