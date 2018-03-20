#pragma once
#ifndef __PIC2D_H__
#define __PIC2D_H__

#include <Eigen/Dense>
#include <vector>
#include "Array2d.h"
#include "FluidWorld2D.h"

using namespace PICArray2d;
using namespace Eigen;
using namespace std;

enum Gid { B = -1, A = 0, F = 1 };

class PIC2D
{
public:
	PIC2D() {};

	void clean();

	void Initialize(FluidWorld2D* p_world, Vector2f p_origin, Vector2f p_bSize, Vector2i p_nCount, float p_rho);
	void AssignBoundary(std::vector<FParticle2D*>& p_BpList);
	void AssignCells(FluidWorld2D* p_world);
	void GetNeigboringParticles_cell(int i, int j, int wl, int wh, int hl, int hh, std::vector<FParticle2D *>& res);
	void Map_P2G(FluidWorld2D* p_world);
	Vector2i GetNiNj() { return Vector2i(ni, nj); }
	Vector2f GetDxDy() { return Vector2f(dx, dy); }

	Vector2f GetVelocity(int p_gridIdx);
	void UpdateAffineMatrix(FluidWorld2D* p_world);

	Vector2i& GetAssignResultF(int p_idx) { return AssignResultF[p_idx]; }
	void GetDescriptorAll(float result[], int p_bound_cnt);
	void GetDescriptor(Vector2i p_ij, float result[], int desc_width);

	int GetGridSize() { return ni*nj; }
	Gid GetGid(int p_idx) { return geo[p_idx]; }
	Vector2f& GetGridPos(int p_idx) { return cells_centor_pos[p_idx]; }
	Vector2f& GetGridPos(int i, int j, int k) { return cells_centor_pos[j*ni + i]; }

private:
	float rho;
	Vector2f origin;
	int ni, nj;
	float dx, dy;

	std::vector<Gid> geo;
	PICArray2d::Array2d<float> u;
	PICArray2d::Array2d<float> v;

	std::vector<std::vector<FParticle2D*>> cellsForF;
	std::vector<std::vector<FParticle2D*>> cellsForB;
	std::vector<Vector2f> cells_centor_pos;
	std::vector<Vector2f> cells_centor_uv_coord;
	std::vector<Vector2i> AssignResultF;

	inline float interpolate_value(Vector2f& point, PICArray2d::Array2d<float>& grid);
	inline Vector2f affine_interpolate_value(Vector2f& point, PICArray2d::Array2d<float>& grid);
};

inline float linear_kernel(const Vector2f& d, const float& h)
{
	return std::max((1.0 - fabs(d(0) / h)) * (1.0 - fabs(d(1) / h)) * (1.0 - fabs(d(2) / h)), 0.0);
}

template<class T>
inline void get_barycentric(T x, int& i, T& f, int i_low, int i_high)
{
	T s = std::floor(x);
	i = (int)s;
	if (i<i_low) {
		i = i_low;
		f = 0;
	}
	else if (i>i_high - 2) {
		i = i_high - 2;
		f = 1;
	}
	else
		f = (T)(x - s);
}

template<class S, class T>
inline S lerp(const S& value0, const S& value1, T f)
{
	return (1 - f)*value0 + f*value1;
}

template<class S, class T>
inline S bilerp(const S& v00, const S& v10,
	const S& v01, const S& v11,
	T fx, T fy)
{
	return lerp(lerp(v00, v10, fx),
		lerp(v01, v11, fx),
		fy);
}

template<class S, class T>
inline S trilerp(const S& v000, const S& v100,
	const S& v010, const S& v110,
	const S& v001, const S& v101,
	const S& v011, const S& v111,
	T fx, T fy, T fz)
{
	return lerp(bilerp(v000, v100, v010, v110, fx, fy),
		bilerp(v001, v101, v011, v111, fx, fy),
		fz);
}

template<class T>
inline Eigen::Matrix<T, 2, 1> grad_bilerp(const T& v00, const T& v10,
	const T& v01, const T& v11,
	T fx, T fy)
{
	return Eigen::Matrix<T, 2, 1>(fy - 1.0, fx - 1.0) * v00 +
		Eigen::Matrix<T, 2, 1>(1.0 - fy, -fx) * v10 +
		Eigen::Matrix<T, 2, 1>(-fy, 1.0 - fx) * v01 +
		Eigen::Matrix<T, 2, 1>(fy, fx) * v11;
}

template<class T>
inline Eigen::Matrix<T, 3, 1> grad_trilerp(const T& v000, const T& v100, const T& v010, const T& v110,
	const T& v001, const T& v101, const T& v011, const T& v111,
	T fx, T fy, T fz)
{

	return  Eigen::Matrix<T, 3, 1>(fy - 1.0, fz - 1.0, fx - 1.0) * v000 +
		Eigen::Matrix<T, 3, 1>(1.0 - fy, fz - 1.0, -fx) * v100 +
		Eigen::Matrix<T, 3, 1>(-fy, 1.0 - fx, fx - 1.0) * v010 +
		Eigen::Matrix<T, 3, 1>(fy, 1.0 - fz, -fx) * v110 +
		Eigen::Matrix<T, 3, 1>(fy - 1.0, -fz, 1.0 - fx) * v001 +
		Eigen::Matrix<T, 3, 1>(1.0 - fy, -fz, fx) * v101 +
		Eigen::Matrix<T, 3, 1>(-fy, fz, 1.0 - fx) * v011 +
		Eigen::Matrix<T, 3, 1>(fy, fz, fx) * v111;
}


#endif