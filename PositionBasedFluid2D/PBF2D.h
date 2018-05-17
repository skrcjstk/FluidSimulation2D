#pragma once
#ifndef __PBF2D_H__
#define __PBF2D_H__

#include <stdlib.h>
#include <Eigen/Dense>
#include <vector>
#include "FParticle2D.h"
#include "FluidKernel2D.h"

using namespace Eigen;
using namespace std;

class PBFWorld2D
{
private:
	// simulation data
	std::vector<float> m_particlesLambda;
	std::vector<Vector2f> m_deltaX;

	float m_viscosity;
	float m_surfaceTensionThr;
	float m_surfaceTensionCoeff;
	float m_restDensity;
	float m_init_sl;

public:
	std::vector<Vector2f> m_deltaFromF;
	std::vector<Vector2f> m_deltaFromB;
	
	PBFWorld2D(float p_restDensity, float p_viscosity, float p_surfaceTensionThr, float p_surfaceTensionCoeff);
	void Reset();

	void InitializeSimulationData(int p_numOfParticles);

	void ComputeXSPHViscosity(std::vector<FParticle2D*>& m_particles, FluidKernel2D& p_kernel);
	void ConstraintProjection(std::vector<FParticle2D*>& p_particles, std::vector<FParticle2D*>& p_boundaryParticles, float p_timeStep, FluidKernel2D& p_kernel);

};

#endif