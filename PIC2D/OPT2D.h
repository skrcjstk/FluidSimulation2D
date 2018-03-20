#pragma once
#ifndef __OPT2D_H__
#define __OPT2D_H__

#include <iostream>
#include <Eigen/dense>
#include <vector>
#include "FluidWorld2D.h"

using namespace Eigen;

class GridCell
{
public:
	GridCell()
	{
		m_centerPosition[0] = m_centerPosition[1] = m_centerPosition[2] = 0.0f;
		m_inFluid = false;
		m_inBoundary = false;
	}

	void SetValue(Vector2f p_cp, float p_gridSize)
	{
		m_centerPosition[0] = p_cp[0];
		m_centerPosition[1] = p_cp[1];

		float halfGridSize = 0.5f * p_gridSize;
		m_minCorner[0] = m_centerPosition[0] - halfGridSize;
		m_minCorner[1] = m_centerPosition[1] - halfGridSize;
		m_maxCorner[0] = m_centerPosition[0] + halfGridSize;
		m_maxCorner[1] = m_centerPosition[1] + halfGridSize;
	}

	Vector2f m_centerPosition;
	Vector2f m_minCorner;
	Vector2f m_maxCorner;

	bool m_inFluid;
	bool m_inBoundary;

	std::vector<Vector2f> m_boundaryNeighborList;

};

class FlowBoundary
{
public:
	FlowBoundary()
	{
		m_gridWidth = 0.0f;
		m_cellCnt = 0;
		m_fThr = 0.0f;
		m_lcp = 0.0f;

	}
	~FlowBoundary()
	{
		m_boundary.clear();
		m_boundary.resize(0);
	}

	void FlowBoundary::SetGridWidth(float p_gridWidth)
	{
		m_gridWidth = p_gridWidth;
	}

	void FlowBoundary::SetFluidThreshold(float p_thr)
	{
		m_fThr = p_thr;
	}

	GridCell& GetGridCell(int p_idx)
	{
		return m_boundary[p_idx];
	}

	void CreateFlowBoundary(Vector2f start, Vector2f end, float p_gridWidth);
	void CheckInBoundary(std::vector<Vector2f>& p_boundaryParticles);
	
	std::vector<GridCell>& GetFlowBoundary() { return m_boundary; }
	int GetCellCnt() { return m_cellCnt; }
	float GetGridWidth() { return m_gridWidth; }


	void CreateCoarsePS(FluidWorld2D* p_mainWorld, std::vector<Vector2f>& p_subWorldFluidParticle);
	void InitializeDataStructure(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld);
	void NeighborSearchBTWTwoRes(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld);

	void InterpolateVelocity(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld);
	void NeighborSearchBTWTwoRes2(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld);

private:
	void CreateBoundaryWall(Vector2f p_min, Vector2f p_max);

	float m_gridWidth;
	float m_lcp;
	int m_cntOfcandisPerOneAxis;

	std::vector<GridCell> m_boundary;
	int m_cellCnt;
	float m_fThr;
	FluidKernel2D k;

	std::vector<std::vector<int>> m_neighborListforSub;
};


#endif