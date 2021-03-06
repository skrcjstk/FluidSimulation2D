#pragma once
#ifndef __FLUID_WORLD2D_H__
#define __FLUID_WORLD2D_H__

#include <vector>
#include <chrono>
#include <Eigen/Dense>
#include "FParticle2D.h"
#include "FluidKernel2D.h"
#include "PBF2D.h"

using namespace Eigen;

class FluidWorld2D
{
private:
	PBFWorld2D* pbfWorld;

	std::vector<FParticle2D*> m_particles;
	std::vector<FParticle2D*> m_boundaryParticles;

	FluidKernel2D m_kernel;

	int  m_numOfParticles;
	int  m_numOfBoundaryParticles;

	float m_fluidParticleRadius;
	float m_restDensity;
	float m_timeStep;
	bool  m_useGravity;

	void NeighborListUpdate();
	void UpdateTimeStepSizeCFL();

public:
	float m_accTimeIntegration;

	FluidWorld2D();
	~FluidWorld2D();

	void Reset();
	void CreateBoundaryParticles(std::vector<Vector2f>& p_particles, float p_particleRadius);
	void CreateFluidParticles(std::vector<Vector2f>& p_particles, float p_particleRadius);

	void StepPBF();
	//void StepPBFonMain1();
	//void StepPBFonMain2();
	void StepPBFonSub1();
	void StepPBFonSub2();
	void StepPBFonSub1WithTF();
	void StepPBFonSub2WithTF();

	PBFWorld2D* GetPBFWorld2D() { return pbfWorld; }

	FParticle2D* GetParticle(int p_index)
	{
		return m_particles[p_index];
	}
	std::vector<FParticle2D*>& GetParticleList()
	{
		return m_particles;
	}
	FParticle2D* GetBoundaryParticle(int p_index)
	{
		return m_boundaryParticles[p_index];
	}
	std::vector<FParticle2D*>& GetBoundaryParticleList()
	{
		return m_boundaryParticles;
	}

	int  GetNumOfParticles()
	{
		return m_numOfParticles;
	}
	int GetNumOfBoundaryParticles()
	{
		return m_numOfBoundaryParticles;
	}
	void SetUseGravity(bool p_useGravity)
	{
		m_useGravity = p_useGravity;
	}
	void SetTimeStep(float p_timeStep)
	{
		m_timeStep = p_timeStep;
	}
	float GetParticleRadius()
	{
		return m_fluidParticleRadius;
	}

	FluidKernel2D& GetKernel() { return m_kernel; }
	float GetTimeStep() { return m_timeStep; }
	float GetFluidParticleRadius() { return m_fluidParticleRadius; }
	float GetRestDensity() { return m_restDensity; }

};


#endif
