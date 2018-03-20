#include "OPT2D.h"

void FlowBoundary::CreateFlowBoundary(Vector2f start, Vector2f end, float p_gridWidth)
{
	m_gridWidth = p_gridWidth;
	m_lcp = m_gridWidth * 0.01f;
	m_cntOfcandisPerOneAxis = int((m_gridWidth - m_lcp) / m_lcp) + 1;

	CreateBoundaryWall(start, end);
}

void FlowBoundary::CreateBoundaryWall(Vector2f p_min, Vector2f p_max)
{
	Vector2f diff = p_max - p_min;

	int countX = (int)(diff[0] / m_gridWidth) + 1;
	int countY = (int)(diff[1] / m_gridWidth) + 1;

	printf("Flow Boundary Wall - diff(%f, %f)\n", diff[0], diff[1]);
	printf("Flow Boundary Wall - count(%d, %d)\n", countX, countY);

	m_boundary.resize(countX*countY);
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		
		for (int j = 0; j < countY; j++)
		{
			for (int i = 0; i < countX; i++)
			{
				Vector2f position = p_min + Vector2f(i*m_gridWidth, j*m_gridWidth);
				m_boundary[j*countX + i].SetValue(position, m_gridWidth);
			}
		}
	}
	m_cellCnt = (int)m_boundary.size();
}

void FlowBoundary::CheckInBoundary(std::vector<Vector2f>& p_boundaryParticles)
{
	float coarseR = m_gridWidth * 0.5f;
	float cellArea = m_gridWidth * m_gridWidth;

	for (int i = 0; i < m_cellCnt; i++) {
		float BoundaryRatio = 0.0f;
		float cell_minX = m_boundary[i].m_minCorner[0];
		float cell_minY = m_boundary[i].m_minCorner[1];
		float cell_maxX = m_boundary[i].m_maxCorner[0];
		float cell_maxY = m_boundary[i].m_maxCorner[1];

		for (int j = 0; j < p_boundaryParticles.size(); j++)
		{
			float p_minX = p_boundaryParticles[j][0] - coarseR;
			float p_minY = p_boundaryParticles[j][1] - coarseR;
			float p_maxX = p_boundaryParticles[j][0] + coarseR;
			float p_maxY = p_boundaryParticles[j][1] + coarseR;

			float x_overlap = std::max(0.0f, std::min(cell_maxX, p_maxX) - std::max(cell_minX, p_minX));
			float y_overlap = std::max(0.0f, std::min(cell_maxY, p_maxY) - std::max(cell_minY, p_minY));
			float overlap_size = x_overlap * y_overlap;

			if (overlap_size > 0.0f)
			{
				BoundaryRatio += overlap_size;
			}
		}

		BoundaryRatio = BoundaryRatio / cellArea;
		if (BoundaryRatio > m_fThr)
		{
			m_boundary[i].m_inBoundary = true;
		}
	}
}

void FlowBoundary::CreateCoarsePS(FluidWorld2D* p_mainWorld, std::vector<Vector2f>& p_subWorldFluidParticle)
{
	std::vector<Vector2f> neighborFineParticles;
	float coarseR = m_gridWidth * 0.5f;
	float fineR = p_mainWorld->GetParticleRadius();
	std::vector<FParticle2D*>& fineP = p_mainWorld->GetParticleList();
	float cellArea = m_gridWidth * m_gridWidth;

	for (int i = 0; i < m_cellCnt; i++) {
		float fluidRatio = 0.0f;
		float cell_minX = m_boundary[i].m_minCorner[0];
		float cell_minY = m_boundary[i].m_minCorner[1];
		float cell_maxX = m_boundary[i].m_maxCorner[0];
		float cell_maxY = m_boundary[i].m_maxCorner[1];

		for (int j = 0; j < fineP.size(); j++)
		{
			float p_minX = fineP[j]->m_curPosition[0] - fineR;
			float p_minY = fineP[j]->m_curPosition[1] - fineR;
			float p_maxX = fineP[j]->m_curPosition[0] + fineR;
			float p_maxY = fineP[j]->m_curPosition[1] + fineR;

			float x_overlap = std::max(0.0f, std::min(cell_maxX, p_maxX) - std::max(cell_minX, p_minX));
			float y_overlap = std::max(0.0f, std::min(cell_maxY, p_maxY) - std::max(cell_minY, p_minY));
			float overlap_size = x_overlap * y_overlap;

			if (overlap_size > 0.0f)
			{
				fluidRatio += overlap_size;
			}
		}

		fluidRatio = fluidRatio / cellArea;
		if (fluidRatio > m_fThr && m_boundary[i].m_inBoundary != true)
		{
			m_boundary[i].m_inFluid = true;
			p_subWorldFluidParticle.push_back(Vector2f(m_boundary[i].m_centerPosition[0], m_boundary[i].m_centerPosition[1]));
		}
	}

	m_neighborListforSub.resize((int)p_subWorldFluidParticle.size());
	
}

void FlowBoundary::NeighborSearchBTWTwoRes(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld)
{
	// weight = w(r,h) = w(r, L0) 
	k.SetSmoothingRadius(m_gridWidth);

	std::vector<FParticle2D*>& mainP = p_mainWorld->GetParticleList();
	std::vector<FParticle2D*>& subP = p_subWorld->GetParticleList();

	for (int i = 0; i < (int)subP.size(); i++)
	{
		m_neighborListforSub[i].clear();
		m_neighborListforSub[i].resize(0);
		Vector2f& subPos = subP[i]->m_curPosition;

		for (int j = 0; j < (int)mainP.size(); j++)
		{
			Vector2f& mainPos = mainP[j]->m_curPosition;

			if ((subPos - mainPos).norm() <= m_gridWidth)
				m_neighborListforSub[i].push_back(j);
		}
	}
}

void FlowBoundary::InterpolateVelocity(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld)
{
	std::vector<FParticle2D*>& mainP = p_mainWorld->GetParticleList();
	std::vector<FParticle2D*>& subP = p_subWorld->GetParticleList();

	// tempVel for coarseP
	float h = p_subWorld->GetTimeStep();

	// MLS interpolation
	VectorXf acc_distVelX(6), acc_distVelY(6);
	VectorXf res_distVelX(6), res_distVelY(6);
	VectorXf distFeature(6);
	MatrixXf acc_distMat(6, 6);

	for (int i = 0; i < (int)subP.size(); i++)
	{
		acc_distMat.setZero();
		acc_distVelX.setZero();
		acc_distVelY.setZero();

		if ((int)m_neighborListforSub[i].size() > 0)
		{
			Vector2f& subPos = subP[i]->m_curPosition;

			for (int j = 0; j < (int)m_neighborListforSub[i].size(); j++)
			{
				FParticle2D* mainNeighbor = mainP[m_neighborListforSub[i][j]];
				Vector2f& mainPos = mainNeighbor->m_curPosition;
				Vector2f dist = mainPos - subPos;
				float weight = k.Cubic_Kernel(dist);

				distFeature[0] = 1;
				distFeature[1] = dist[0];
				distFeature[2] = dist[1];
				distFeature[3] = dist[0] * dist[0];
				distFeature[4] = dist[0] * dist[1];
				distFeature[5] = dist[1] * dist[1];

				for (int m = 0; m < 6; m++)
					for (int n = 0; n < 6; n++)
						acc_distMat(m, n) += weight * distFeature[m] * distFeature[n];

				for (int m = 0; m < 6; m++)
				{
					acc_distVelX[m] += weight * distFeature[m] * mainNeighbor->m_velocity[0];
					acc_distVelY[m] += weight * distFeature[m] * mainNeighbor->m_velocity[1];
				}
			}

			// velocity interpolation
			res_distVelX = acc_distMat.jacobiSvd(ComputeThinU | ComputeThinV).solve(acc_distVelX);
			res_distVelY = acc_distMat.jacobiSvd(ComputeThinU | ComputeThinV).solve(acc_distVelY);
			//printf("(%d) res_distVel(%f, %f, %f)\n", i, res_distVelX[0], res_distVelY[0], res_distVelZ[0]);

			subP[i]->m_velocity[0] = res_distVelX[0];
			subP[i]->m_velocity[1] = res_distVelY[0];
		}
		else
			subP[i]->m_velocity = Vector2f(0.0f, 0.0f);

	}
}