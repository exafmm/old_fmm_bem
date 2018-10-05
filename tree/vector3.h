/*
 *  vector3.h: 3D vector3 operations include file
 *.............................................................................
 *    version 1:  Dec 1992   Piet Hut, Steve McMillan, Jun Makino
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of class vector3
 *.............................................................................
 */

// This is slightly modified version of vector3 header file
// taken from STARLAB
// -- J. Makino

#ifndef STARLAB_vector3_H
#define STARLAB_vector3_H

#include <iostream>
#include <cmath>

inline real readjust_r(real x)
{
  return x;
}

/*-----------------------------------------------------------------------------
 *  vector3  --  a class for 3-dimensional vector3s
 *-----------------------------------------------------------------------------
 */

const int ndim = 3;

class vector3
{
  private:

    real x, y, z;

  public:

    //  Default: initialize to zero.
    inline vector3(real c = 0) : x(c), y(c), z(c) {}
    inline vector3(real _x, real _y, real _z) : x(_x), y(_y), z(_z) {}
    //  []: the return type is declared as a reference (&), so that it can be used
    //  on the left-hand side of an asignment, as well as on the right-hand side,
    //  i.e.  v[1] = 3.14  and  x = v[2]  are both allowed and work as expected.
    inline real & operator [] (int i)       {return (&x)[i];}
    inline const real & operator [] (int i) const { return (&x)[i]; }
    inline void print() {cout << x << " " << y << " " << z << "\n";}
    inline vector3  readjust(){
      return vector3(*this);
    }

    bool has_nan() const {
      return __builtin_isnan(x)
        || __builtin_isnan(y)
        || __builtin_isnan(z);
    }
    //  Unary -
    inline vector3 operator - ()
    { return vector3(-x, -y, -z); }
    //  Dot product.
    inline real operator * (const vector3& b) const
    { return x * b.x + y * b.y + z * b.z; }
    //  Outer product.
    inline vector3 operator ^ (const vector3 &b)
    { return vector3(y*b.z - z*b.y,
                     z*b.x - x*b.z,
                     x*b.y - y*b.x);}
    //  vector3 +, -
    inline vector3 operator + (const vector3 &b)
    { return vector3(x+b.x, y+b.y, z+b.z); }
    inline vector3 operator - (const vector3 &b)
    { return vector3(x-b.x, y-b.y, z-b.z); }
    friend vector3 operator + (real, const vector3 & );
    friend vector3 operator + (const vector3 &, real);
    //  Scalar *, /
    friend vector3 operator * (real, const vector3 & );
    friend vector3 operator * (const vector3 &, real);
    friend vector3 operator / (const vector3 &, real);
    //  vector3 +=, -=, *=, /=
    inline vector3& operator += (const vector3& b)
    {
      x += b.x;       
      y += b.y;
      z += b.z;
      return *this;
    }
    inline vector3& operator -= (const vector3& b)
    {
      x -= b.x;
      y -= b.y;
      z -= b.z;
      return *this;
    }
    inline vector3& operator *= (const real b)
    {
      x *= b;
      y *= b;
      z *= b;
      return *this;
    }
    inline vector3& operator /= (const real b)
    {
      register real binv = 1.0/b;
      x *= binv;
      y *= binv;
      z *= binv;
      return *this;
    }
    //      Input / Output
    friend ostream & operator << (ostream & , const vector3 & );
    friend istream & operator >> (istream & , vector3 & );
    real square() const { return (*this)*(*this); }
    real abs()    const { return std::sqrt(square()); }
};

inline  ostream & operator << (ostream & s, const vector3 & v)
{ return s << v.x << "  " << v.y << "  " << v.z; }

inline  istream & operator >> (istream & s, vector3 & v)
{
  s >> v.x >> v.y >> v.z;
  return s;
}

inline  real square(const vector3 v) {return v*v;}
inline  real abs   (const vector3 v) {return sqrt(v*v);}

inline  vector3 operator + (real b, const vector3 & v)
{ return vector3(b+v.x, b+v.y, b+v.z); }

inline  vector3 operator + (const vector3 & v, real b)
{ return vector3(b+v.x, b+v.y, b+v.z); }

inline  vector3 operator * (real b, const vector3 & v)
{ return vector3(b*v.x, b*v.y, b*v.z); }

inline  vector3 operator * (const vector3 & v, real b)
{ return vector3(b*v.x, b*v.y, b*v.z); }

inline  vector3 operator / (const vector3 & v, real b)
{ return vector3(v.x/b, v.y/b, v.z/b); }

#endif

