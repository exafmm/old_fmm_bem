// Copyright(C) 2010 by 
// Tsuyoshi Hamada <hamada@progrape.jp>
// Keigo Nitadori <nitadori@margaux.astron.s.u-tokyo.ac.jp>
// Rio Yokota <rio.yokota@bristol.ac.uk>

#ifndef __V4SF_H
#define __V4SF_H
struct v4sf{
	typedef float _v4sf __attribute__ ((vector_size(16)));
	_v4sf val;

	v4sf() {}
	v4sf(float f) : val((_v4sf){f, f, f, f}) {}
	v4sf(_v4sf v): val(v) {}

	const v4sf abs() const{
		typedef int v4si __attribute__ ((vector_size(16)));
		v4si mask = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};
		return __builtin_ia32_andps(val, (_v4sf)mask);
	}

	const v4sf operator+(const v4sf &rhs) const{
		return val + rhs.val;
	}
	const v4sf operator-(const v4sf &rhs) const{
		return val - rhs.val;
	}
	const v4sf operator*(const v4sf &rhs) const{
		return val * rhs.val;
	}
	const v4sf operator/(const v4sf &rhs) const{
		return val / rhs.val;
	}
	const v4sf operator+=(const v4sf &rhs){
		val += rhs.val;
		return val;
	}
	const v4sf operator-=(const v4sf &rhs){
		val -= rhs.val;
		return val;
	}
	const v4sf operator*=(const v4sf &rhs){
		val *= rhs.val;
		return val;
	}
	const v4sf operator/=(const v4sf &rhs){
		val /= rhs.val;
		return val;
	}
	const v4sf operator<(const v4sf &rhs) const{
		return (_v4sf)__builtin_ia32_cmpltps(val, rhs.val);
	}
	const v4sf operator<=(const v4sf &rhs) const{
		return (_v4sf)__builtin_ia32_cmpleps(val, rhs.val);
	}
	const v4sf operator>(const v4sf &rhs) const{
		return (_v4sf)__builtin_ia32_cmpgtps(val, rhs.val);
	}
	const v4sf operator>=(const v4sf &rhs) const{
		return (_v4sf)__builtin_ia32_cmpgeps(val, rhs.val);
	}
	const v4sf operator|(const v4sf &rhs) const{
		return __builtin_ia32_orps(val, rhs.val);
	}
	const v4sf operator&(const v4sf &rhs) const{
		return __builtin_ia32_andps(val, rhs.val);
	}
};

struct v4sf_x4{
	typedef float _v4sf __attribute__ ((vector_size(16)));
	_v4sf x, y, z, w;
	v4sf_x4() {}
	v4sf_x4(_v4sf v0, _v4sf v1, _v4sf v2, _v4sf v3, bool transpose = false){
		if(!transpose){
			x = v0;
			y = v1;
			z = v2;
			w = v3;
		}else{
			_v4sf t0 = __builtin_ia32_unpcklps(v0, v2);
			_v4sf t1 = __builtin_ia32_unpckhps(v0, v2);
			_v4sf t2 = __builtin_ia32_unpcklps(v1, v3);
			_v4sf t3 = __builtin_ia32_unpckhps(v1, v3);
			x = __builtin_ia32_unpcklps(t0, t2);
			y = __builtin_ia32_unpckhps(t0, t2);
			z = __builtin_ia32_unpcklps(t1, t3);
			w = __builtin_ia32_unpckhps(t1, t3);
		}
	}
	v4sf_x4(_v4sf v){
		x = __builtin_ia32_shufps(v, v, 0x00);
		y = __builtin_ia32_shufps(v, v, 0x55);
		z = __builtin_ia32_shufps(v, v, 0xaa);
		w = __builtin_ia32_shufps(v, v, 0xff);
	}
	v4sf_x4 transpose(){
		return v4sf_x4(x, y, z, w, true);
	}
};
#endif // __V4SF_H
