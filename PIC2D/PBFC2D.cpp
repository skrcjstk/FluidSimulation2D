#include "PBFC2D.h"
PBFC2D::PBFC2D(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld)
{
	float mainR = p_mainWorld->GetParticleRadius();
	constKernel.SetSmoothingRadius(8.0f * mainR);

	m_PBFCData.m_lambdaForMain.resize(p_mainWorld->GetNumOfParticles());
	m_PBFCData.m_corrWithDensity.resize(p_subWorld->GetNumOfParticles());
	m_PBFCData.m_corrWithVelocity.resize(p_subWorld->GetNumOfParticles());
	m_PBFCData.m_weightForVelocityC.resize(p_subWorld->GetNumOfParticles());

	m_neighListwithSubP.resize(p_mainWorld->GetNumOfParticles());
	m_neighListwithSubBoundaryP.resize(p_mainWorld->GetNumOfParticles());

	m_tDataForCoarse.resize(p_mainWorld->GetNumOfParticles());
	m_tDataForFine.resize(p_mainWorld->GetNumOfParticles());
	m_tDataForCoarseB.resize(p_mainWorld->GetNumOfParticles());
	m_tDataForFineB.resize(p_mainWorld->GetNumOfParticles());
}
void PBFC2D::NeighborBTWTwoResForPBFC(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld)
{
	float searchRange = constKernel.GetSmoothingRadius();
	std::vector<FParticle2D*>& mainP = p_mainWorld->GetParticleList();
	std::vector<FParticle2D*>& subP = p_subWorld->GetParticleList();
	std::vector<FParticle2D*>& boundaryMainP = p_mainWorld->GetBoundaryParticleList();
	std::vector<FParticle2D*>& boundarySubP = p_subWorld->GetBoundaryParticleList();
	
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < mainP.size(); i++)
		{
			Vector2f& mainPos = mainP[i]->m_curPosition;

			m_neighListwithSubP[i].clear();
			m_neighListwithSubP[i].resize(0);
			for (int j = 0; j < subP.size(); j++)
			{
				Vector2f& subPos = subP[j]->m_curPosition;
				if ((mainPos - subPos).norm() <= searchRange)
					m_neighListwithSubP[i].push_back(j);
			}
			m_neighListwithSubBoundaryP[i].clear();
			m_neighListwithSubBoundaryP[i].resize(0);
			for (int j = 0; j < boundarySubP.size(); j++)
			{
				Vector2f& subPos = boundarySubP[j]->m_curPosition;
				if ((mainPos - subPos).norm() <= searchRange)
					m_neighListwithSubBoundaryP[i].push_back(j);
			}
		}
	}
}
void PBFC2D::SolvePBFCConstaints(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld)
{
	int numOfMainP = p_mainWorld->GetNumOfParticles();
	int numOfSubP = p_subWorld->GetNumOfParticles();
	std::vector<FParticle2D*>& mainP = p_mainWorld->GetParticleList();
	std::vector<FParticle2D*>& subP = p_subWorld->GetParticleList();
	std::vector<FParticle2D*>& boundarySubP = p_subWorld->GetBoundaryParticleList();

	float mainRestDensity = p_mainWorld->GetRestDensity();

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < numOfSubP; i++)
		{
			m_PBFCData.m_corrWithDensity[i].setZero();
			m_PBFCData.m_corrWithVelocity[i].setZero();
			m_PBFCData.m_weightForVelocityC[i] = 0.0f;
		}

		// update CoarseLambda & correction with Density Constraint
#pragma omp for schedule(static) 
		for (int i = 0; i < numOfMainP; i++)
		{
			m_PBFCData.m_lambdaForMain[i] = 0.0f;
			Vector2f& mainPos = mainP[i]->m_curPosition;

			float density = 0.0f;
			for (int j = 0; j < m_neighListwithSubP[i].size(); j++)
			{
				int idx = m_neighListwithSubP[i][j];
				Vector2f& subPos = subP[idx]->m_curPosition;
				density += subP[idx]->m_mass * constKernel.Cubic_Kernel(mainPos - subPos);
			}
			for (int j = 0; j < m_neighListwithSubBoundaryP[i].size(); j++)
			{
				int idx = m_neighListwithSubBoundaryP[i][j];
				Vector2f& subPos = boundarySubP[idx]->m_curPosition;
				density += boundarySubP[idx]->m_mass * constKernel.Cubic_Kernel(mainPos - subPos);
			}

			float C = std::max(density / mainRestDensity - 1.0f, 0.0f);

			if (C != 0.0f)
			{
				// Compute gradients dC/dx_j 
				float sum_grad_C2 = 0.0;
				Vector2f gradC_i(0.0f, 0.0f);

				for (int j = 0; j < m_neighListwithSubP[i].size(); j++)
				{
					int idx = m_neighListwithSubP[i][j];
					Vector2f& subPos = subP[idx]->m_curPosition;

					Vector2f gradC_j = -subP[idx]->m_mass / mainRestDensity * constKernel.Cubic_Kernel_Gradient(mainPos - subPos);
					sum_grad_C2 += gradC_j.squaredNorm();
					gradC_i -= gradC_j;
				}

				for (int j = 0; j < m_neighListwithSubBoundaryP[i].size(); j++)
				{
					int idx = m_neighListwithSubBoundaryP[i][j];
					Vector2f& subPos = boundarySubP[idx]->m_curPosition;

					Vector2f gradC_j = -boundarySubP[idx]->m_mass / mainRestDensity * constKernel.Cubic_Kernel_Gradient(mainPos - subPos);
					sum_grad_C2 += gradC_j.squaredNorm();
					gradC_i -= gradC_j;
				}

				sum_grad_C2 += gradC_i.squaredNorm();

				// Compute lambda
				m_PBFCData.m_lambdaForMain[i] = -C / (sum_grad_C2 + 1.0e-6);
			}

			// calc correction with density constraint
			if (m_PBFCData.m_lambdaForMain[i] != 0.0f)
			{
				for (int j = 0; j < m_neighListwithSubP[i].size(); j++)
				{
					int idx = m_neighListwithSubP[i][j];
					Vector2f& subPos = subP[idx]->m_curPosition;

					Vector2f gradC_j = -subP[idx]->m_mass / mainRestDensity * constKernel.Cubic_Kernel_Gradient(mainPos - subPos);
					m_PBFCData.m_corrWithDensity[idx] += m_intensityOfDensityC * m_PBFCData.m_lambdaForMain[i] * gradC_j;
				}
			}

			// calc correction with Velocity Constraint part1
			for (int j = 0; j < m_neighListwithSubP[i].size(); j++)
			{
				int idx = m_neighListwithSubP[i][j];
				Vector2f& subPos = subP[idx]->m_curPosition;

				m_PBFCData.m_corrWithVelocity[idx] += mainP[i]->m_velocity * constKernel.Cubic_Kernel(mainPos - subPos);
				m_PBFCData.m_weightForVelocityC[idx] += constKernel.Cubic_Kernel(mainPos - subPos);
			}
		}

#pragma omp for schedule(static) 
		// calc correction with Velocity Constraint part2
		for (int i = 0; i < numOfSubP; i++)
		{
			if (m_PBFCData.m_corrWithVelocity[i].norm() != 0.0f)
			{
				m_PBFCData.m_corrWithVelocity[i] = m_PBFCData.m_corrWithVelocity[i] / m_PBFCData.m_weightForVelocityC[i];
				m_PBFCData.m_corrWithVelocity[i] = m_intensityOfVelocityC * p_subWorld->GetTimeStep() * (m_PBFCData.m_corrWithVelocity[i] - subP[i]->m_velocity);
			}
			subP[i]->m_curPosition += m_PBFCData.m_corrWithDensity[i] + m_PBFCData.m_corrWithVelocity[i];
			//m_deltaPWithControl[i] = m_PBFCData.m_corrWithVelocity[i] + m_PBFCData.m_corrWithDensity[i];
		}
	}
}

void PBFC2D::UpdateTrainingData(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld)
{
	float fineR = p_mainWorld->GetParticleRadius();
	float searchRange = 8.0f * fineR;
	float halfRange = 4.0f * fineR;

	std::vector<FParticle2D*>& fineP = p_mainWorld->GetParticleList();
	std::vector<FParticle2D*>& fineBP = p_mainWorld->GetBoundaryParticleList();
	std::vector<FParticle2D*>& coarseP = p_subWorld->GetParticleList();
	std::vector<FParticle2D*>& coarseBP = p_subWorld->GetBoundaryParticleList();
	
	for (int i = 0; i < fineP.size(); i++)
	{
		m_tDataForCoarse[i].clear();
		m_tDataForCoarse[i].resize(0);
		//m_tDataForCoarseB[i].clear();
		//m_tDataForCoarseB[i].resize(0);
		m_tDataForFine[i].clear();
		m_tDataForFine[i].resize(0);
		//m_tDataForFineB[i].clear();
		//m_tDataForFineB[i].resize(0);
	}

	for (int i = 0; i < fineP.size(); i++)
	{
		Vector2f& finePos = fineP[i]->m_oldPosition;
		
		for (int j = 0; j < coarseP.size(); j++)
		{
			Vector2f& coarsePos = coarseP[j]->m_oldPosition;
			Vector2f r = coarsePos - finePos;
			if (r.norm() <= searchRange)
			{
				TrainData a;
				a.mass = coarseP[j]->m_mass;
				a.RVec = r;
				a.dPos = coarseP[j]->m_velocity;
				m_tDataForCoarse[i].push_back(a);
			}
		}
		for (int j = 0; j < coarseBP.size(); j++)
		{
			Vector2f& coarseBPos = coarseBP[j]->m_restPosition;
			Vector2f r = coarseBPos - finePos;
			if (r.norm() <= searchRange)
			{
				TrainData a;
				a.mass = coarseBP[j]->m_mass;
				a.RVec = r;
				a.dPos[0] = a.dPos[1] = 0.0f;
				m_tDataForCoarse[i].push_back(a);
			}
		}
				
		for (int j = 0; j < fineP.size(); j++)
		{
			Vector2f& finePos2 = fineP[j]->m_oldPosition;
			Vector2f r = finePos2 - finePos;
			if (i != j && r.norm() <= halfRange)
			{
				TrainData a;
				a.mass = fineP[j]->m_mass;
				a.RVec = r;
				a.dPos = fineP[j]->m_tempVelocity;
				m_tDataForFine[i].push_back(a);
			}
		}
		for (int j = 0; j < fineBP.size(); j++)
		{
			Vector2f& fineBPos = fineBP[j]->m_restPosition;
			Vector2f r = fineBPos - finePos;
			if (i != j && r.norm() <= halfRange)
			{
				TrainData a;
				a.mass = fineBP[j]->m_mass;
				a.RVec = r;
				a.dPos[0] = a.dPos[1] = 0.0f;
				m_tDataForFine[i].push_back(a);
			}
		}	
	}
}

void PBFC2D::UpdateTrainingDataOnTF(FluidWorld2D* p_mainWorld, FluidWorld2D* p_subWorld)
{
	float fineR = p_mainWorld->GetParticleRadius();
	float searchRange = 8.0f * fineR;
	float halfRange = 4.0f * fineR;

	std::vector<FParticle2D*>& fineP = p_mainWorld->GetParticleList();
	std::vector<FParticle2D*>& fineBP = p_mainWorld->GetBoundaryParticleList();
	std::vector<FParticle2D*>& coarseP = p_subWorld->GetParticleList();
	std::vector<FParticle2D*>& coarseBP = p_subWorld->GetBoundaryParticleList();
	
	for (int i = 0; i < fineP.size(); i++)
	{
		m_tDataForCoarse[i].clear();
		m_tDataForCoarse[i].resize(0);
		//m_tDataForCoarseB[i].clear();
		//m_tDataForCoarseB[i].resize(0);
		m_tDataForFine[i].clear();
		m_tDataForFine[i].resize(0);
		//m_tDataForFineB[i].clear();
		//m_tDataForFineB[i].resize(0);
	}

	for (int i = 0; i < fineP.size(); i++)
	{
		Vector2f& finePos = fineP[i]->m_curPosition;

		for (int j = 0; j < coarseP.size(); j++)
		{
			Vector2f& coarsePos = coarseP[j]->m_oldPosition;
			Vector2f r = coarsePos - finePos;
			if (r.norm() <= searchRange)
			{
				TrainData a;
				a.mass = coarseP[j]->m_mass;
				a.RVec = r;
				a.dPos = coarseP[j]->m_velocity;
				m_tDataForCoarse[i].push_back(a);
			}
		}
		for (int j = 0; j < coarseBP.size(); j++)
		{
			Vector2f& coarseBPos = coarseBP[j]->m_restPosition;
			Vector2f r = coarseBPos - finePos;
			if (r.norm() <= searchRange)
			{
				TrainData a;
				a.mass = coarseBP[j]->m_mass;
				a.RVec = r;
				a.dPos[0] = a.dPos[1] = 0.0f;
				m_tDataForCoarse[i].push_back(a);
			}
		}

		for (int j = 0; j < fineP.size(); j++)
		{
			Vector2f& finePos2 = fineP[j]->m_curPosition;
			Vector2f r = finePos2 - finePos;
			if (i != j && r.norm() <= halfRange)
			{
				TrainData a;
				a.mass = fineP[j]->m_mass;
				a.RVec = r;
				a.dPos = fineP[j]->m_velocity;
				m_tDataForFine[i].push_back(a);
			}
		}
		for (int j = 0; j < fineBP.size(); j++)
		{
			Vector2f& fineBPos = fineBP[j]->m_restPosition;
			Vector2f r = fineBPos - finePos;
			if (i != j && r.norm() <= halfRange)
			{
				TrainData a;
				a.mass = fineBP[j]->m_mass;
				a.RVec = r;
				a.dPos[0] = a.dPos[1] = 0.0f;
				m_tDataForFine[i].push_back(a);
			}
		}
	}
}
