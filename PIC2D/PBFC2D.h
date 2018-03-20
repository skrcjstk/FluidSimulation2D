#pragma once
#ifndef __PBFC2D_H__
#define __PBFC2D_H__

#include <Eigen/dense>
#include <vector>
#include "FluidWorld2D.h"
using namespace Eigen;

class TrainData
{
public:
	int idx;
	float mass;
	float kWeight;
	Vector2f RVec;
	Vector2f kGrad;
	Vector2f dPos;
};

class PBFControlData
{
public:
	std::vector<float> m_lambdaForMain;
	std::vector<Vector2f> m_corrWithDensity;
	std::vector<float> m_weightForVelocityC;
	std::vector<Vector2f> m_corrWithVelocity;
};

class PBFC2D
{
public:
	PBFC2D() {}
	PBFC2D(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld);
	~PBFC2D() {}

	FluidKernel2D constKernel;

	// for TrainingData
	std::vector<std::vector<TrainData>> m_tDataForFine;
	std::vector<std::vector<TrainData>> m_tDataForFineB;
	std::vector<std::vector<TrainData>> m_tDataForCoarse;	
	std::vector<std::vector<TrainData>> m_tDataForCoarseB;

	// for PBFC
	PBFControlData m_PBFCData;
	std::vector<std::vector<int>> m_neighListwithSubP;
	std::vector<std::vector<int>> m_neighListwithSubBoundaryP;

	void NeighborBTWTwoResForPBFC(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld);
	void SolvePBFCConstaints(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld);
	void UpdateTrainingData(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld);
	void UpdateTrainingDataOnTF(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld);

	//void UpdateTrainingDataForSub(FluidWorld2D* p_subWorld);

	float m_intensityOfDensityC = 1.0f;
	float m_intensityOfVelocityC = 1.0f;

};
#endif 