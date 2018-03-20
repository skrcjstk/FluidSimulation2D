#pragma once
#ifndef __FLUID_PARTICLE2D_H__
#define __FLUID_PARTICLE2D_H__

#include <Eigen/Dense>
#include <vector>
#include <utility>

using namespace Eigen;

enum Pid { Fluid, Boundary };

class FParticle2D
{
public:
	Pid   m_pid;
	int   m_pIdx;

	float m_mass;
	float m_density;

	Vector2f m_restPosition;

	Vector2f m_oldPosition;
	Vector2f m_curPosition;
	Vector2f m_tempPosition;
		
	Vector2f m_oldVelocity;
	Vector2f m_velocity;
	Vector2f m_tempVelocity;

	Vector2f m_acceleration;

	std::vector<FParticle2D *> m_neighborList;


	bool m_interpolated = false;

	Matrix2f m_c;
};

#endif
