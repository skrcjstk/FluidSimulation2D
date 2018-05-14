#include "PBF2D.h"

PBFWorld2D::PBFWorld2D(float p_restDensity, float p_viscosity, float p_surfaceTensionThr, float p_surfaceTensionCoeff)
{
	m_restDensity = p_restDensity;
	m_viscosity = p_viscosity;
	m_surfaceTensionThr = p_surfaceTensionThr;
	m_surfaceTensionCoeff = p_surfaceTensionCoeff;
}

void PBFWorld2D::Reset()
{
	int numOfParticles = (int)m_particlesLambda.size();
	for (int i = 0; i < numOfParticles; i++)
	{
		m_particlesLambda[i] = 0.0f;
		m_deltaX[i].setZero();
	}
}

void PBFWorld2D::InitializeSimulationData(int p_numOfParticles)
{
	m_particlesLambda.clear();
	m_particlesLambda.resize(p_numOfParticles);

	m_deltaX.clear();
	m_deltaX.resize(p_numOfParticles);

	m_deltaFromF.clear();
	m_deltaFromF.resize(p_numOfParticles);
	m_deltaFromB.clear();
	m_deltaFromB.resize(p_numOfParticles);
}

void PBFWorld2D::ConstraintProjection(std::vector<FParticle2D*>& p_particles, std::vector<FParticle2D* >& p_boundaryParticles, float p_timeStep)
{
	int maxiter = 100;
	int iter = 0;

	float eps = 1.0e-6f;

	int numParticles = (int)p_particles.size();

	float invH = 1.0f / p_timeStep;

	float density0 = m_restDensity;
	float maxError = 0.01f;
	float eta = maxError * 0.01f * density0;  // maxError is given in percent

	float avg_density_err = 0.0f;

	for (int i = 0; i < numParticles; i++)
	{
		m_deltaFromF[i].setZero();
		m_deltaFromB[i].setZero();
	}

	while (((avg_density_err > eta) || (iter < 2)) && (iter < maxiter))
	{
		avg_density_err = 0.0f;

#pragma omp parallel default(shared)
		{
#pragma omp for schedule(static)
			for (int i = 0; i < numParticles; i++)
			{
				// computePBFDensity
				FParticle2D* pi = p_particles[i];
				pi->m_density = pi->m_mass * m_kernel.Cubic_Kernel0();

				for (unsigned int j = 0; j < pi->m_neighborList.size(); j++)
				{
					FParticle2D* pj = pi->m_neighborList[j];
					Vector2f r = pi->m_curPosition - pj->m_curPosition;

					pi->m_density += pj->m_mass * m_kernel.Cubic_Kernel(r);
				}

				float density_err = std::max(pi->m_density, density0) - density0;
#pragma omp atomic
				avg_density_err += density_err / (float)numParticles;

				// Evaluate constraint function
				float C = std::max(pi->m_density / density0 - 1.0f, 0.0f);

				if (C != 0.0f)
				{
					// Compute gradients dC/dx_j 
					float sum_grad_C2 = 0.0;
					Vector2f gradC_i(0.0f, 0.0f);

					for (unsigned int j = 0; j < pi->m_neighborList.size(); j++)
					{
						FParticle2D* pj = pi->m_neighborList[j];
						Vector2f r = pi->m_curPosition - pj->m_curPosition;

						Vector2f gradC_j = -pj->m_mass / m_restDensity * m_kernel.Cubic_Kernel_Gradient(r);
						sum_grad_C2 += gradC_j.squaredNorm();
						gradC_i -= gradC_j;
					}

					sum_grad_C2 += gradC_i.squaredNorm();

					// Compute lambda
					m_particlesLambda[i] = -C / (sum_grad_C2 + eps);
					//printf("%d's particle C(%f), Lambda(%f)\n", i, C, m_particlesLambda[i]);
				}
				else
				{
					m_particlesLambda[i] = 0.0f;
				}
			}

#pragma omp for schedule(static)
			// Compute position correction
			for (int i = 0; i < numParticles; i++)
			{
				Vector2f corr_fromF(0.0f, 0.0f);
				Vector2f corr_fromB(0.0f, 0.0f);
				FParticle2D* pi = p_particles[i];

				for (unsigned int j = 0; j < pi->m_neighborList.size(); j++)
				{
					FParticle2D* pj = pi->m_neighborList[j];
					Vector2f r = pi->m_curPosition - pj->m_curPosition;
					Vector2f gradC_j = -pj->m_mass / m_restDensity * m_kernel.Cubic_Kernel_Gradient(r);

					if (pj->m_pid == Fluid)
						corr_fromF -= (m_particlesLambda[i] + m_particlesLambda[pj->m_pIdx]) * gradC_j;
					else
						corr_fromB -= (m_particlesLambda[i]) * gradC_j;
				}
				m_deltaX[i] = corr_fromF + corr_fromB;
				
				m_deltaFromF[i] += corr_fromF;
				m_deltaFromB[i] += corr_fromB;
				//if(corr.norm() > 0.0f)
				//	printf("corr(%f, %f)\n", corr[0], corr[1]);
			}

#pragma omp for schedule(static)
			for (int i = 0; i < numParticles; i++)
			{
				p_particles[i]->m_curPosition += m_deltaX[i];
			}
		}

		iter++;
	}
}

void PBFWorld2D::ComputeXSPHViscosity(std::vector<FParticle2D*>& p_particles)
{
	int numParticles = (int)p_particles.size();
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		// Compute viscosity forces (XSPH)
		for (int i = 0; i < numParticles; i++)
		{
			FParticle2D* pi = p_particles[i];
			for (unsigned int j = 0; j < pi->m_neighborList.size(); j++)
			{
				FParticle2D* pj = pi->m_neighborList[j];
				if (pj->m_pid == Fluid)
				{
					Vector2f r = pi->m_curPosition - pj->m_curPosition;
					Vector2f velCorr = m_viscosity * (pj->m_mass / pj->m_density) * (pi->m_velocity - pj->m_velocity) * m_kernel.Cubic_Kernel(r);

					pi->m_velocity -= velCorr;
				}
			}
		}
	}
}