#include "FluidWorld2D.h"

FluidWorld2D::FluidWorld2D()
{
	m_accTimeIntegration = 0.0f;
	m_timeStep = 0.005f;
	m_useGravity = false;
	m_restDensity = 1000.0f;

	pbfWorld = new PBFWorld2D(m_restDensity, 0.02f, 0.2f, 0.2f);
	m_numOfParticles = m_numOfBoundaryParticles = 0;
}
FluidWorld2D::~FluidWorld2D() {}

void FluidWorld2D::CreateBoundaryParticles(std::vector<Vector2f>& p_particles, float p_particleRadius)
{
	FluidKernel2D bKernel;
	float sl = 4.0f * p_particleRadius;
	float diameter = 2.0f * p_particleRadius;
	
	bKernel.SetSmoothingRadius(sl);
	
	int startIndex = (int)m_boundaryParticles.size();
	int nBoundaryParticles = (int)p_particles.size();
	m_boundaryParticles.resize(startIndex + nBoundaryParticles);
	

	// copy boundary particles
#pragma omp for schedule(static)
		for (int i = startIndex; i < (int)m_boundaryParticles.size(); i++)
		{
			m_boundaryParticles[i] = new FParticle2D();
			m_boundaryParticles[i]->m_pid = Boundary;
			m_boundaryParticles[i]->m_pIdx = i;
			m_boundaryParticles[i]->m_restPosition = p_particles[i];
			m_boundaryParticles[i]->m_curPosition = p_particles[i];
			m_boundaryParticles[i]->m_velocity.setZero();
			m_boundaryParticles[i]->m_acceleration.setZero();
		}

	// boudary particles Psi value 
#pragma omp for schedule(static)
	for (int i = startIndex; i < (int)m_boundaryParticles.size(); i++)
	{
		FParticle2D* pi = m_boundaryParticles[i];
		float delta = bKernel.Cubic_Kernel0();
		for (int j = startIndex; j < (int)m_boundaryParticles.size(); j++)
		{
			FParticle2D* pj = m_boundaryParticles[j];
			Vector2f r = pi->m_restPosition - pj->m_restPosition;

			if (i != j && r.norm() <= sl)
				delta += bKernel.Cubic_Kernel(r);
		}
		float volume = 1.0f / delta;
		pi->m_mass = m_restDensity * volume;
	}

	m_numOfBoundaryParticles += nBoundaryParticles;
	printf("Num of boundary particles : %d\n", m_numOfBoundaryParticles);
}

void FluidWorld2D::CreateFluidParticles(std::vector<Vector2f>& p_particles, float p_particleRadius)
{
	int nFluidParticles = (int)p_particles.size();
	m_fluidParticleRadius = p_particleRadius;
	float sl = 4.0f * m_fluidParticleRadius;
	m_kernel.SetSmoothingRadius(sl);
	float diameter = 2.0f * m_fluidParticleRadius;

	int startIndex = (int)m_particles.size();
	m_particles.resize(startIndex + nFluidParticles);
	
#pragma omp parallel default(shared)
	{
		// dam particles creation
#pragma omp for schedule(static)
		for (int i = startIndex; i < (int)m_particles.size(); i++)
		{
			m_particles[i] = new FParticle2D();
			m_particles[i]->m_pid = Fluid;
			m_particles[i]->m_pIdx = i;
			//m_particles[i]->m_mass = 0.8f * m_restDensity * diameter * diameter * diameter; (in 3D case)
			m_particles[i]->m_mass = m_restDensity * p_particleRadius * p_particleRadius;
			m_particles[i]->m_restPosition = p_particles[i];
			m_particles[i]->m_curPosition = p_particles[i];
			m_particles[i]->m_velocity.setZero();
			m_particles[i]->m_acceleration.setZero();
		}
	}
	
	m_numOfParticles += (int)p_particles.size();
	printf("Num of fluid particles : %d\n", m_numOfParticles);

	pbfWorld->InitializeSimulationData(m_numOfParticles);
}

void FluidWorld2D::Reset()
{
	for (int i = 0; i < m_numOfParticles; i++)
	{
		m_particles[i]->m_acceleration.setZero();
		m_particles[i]->m_velocity.setZero();
		m_particles[i]->m_curPosition = m_particles[i]->m_restPosition;
	}

	pbfWorld->Reset();
}
void FluidWorld2D::NeighborListUpdate()
{
	float sl = 4.0f * m_fluidParticleRadius;

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)
		for (int i = 0; i < m_numOfParticles; i++)
		{
			FParticle2D* pi = m_particles[i];
			pi->m_neighborList.clear();
			pi->m_neighborList.resize(0);

			for (int j = 0; j < m_numOfParticles; j++)
			{
				FParticle2D* pj = m_particles[j];
				Vector2f r = pi->m_curPosition - pj->m_curPosition;
				if (i != j && r.norm() < sl)
					pi->m_neighborList.push_back(pj);
			}
			for (int j = 0; j < m_numOfBoundaryParticles; j++)
			{
				FParticle2D* pj = m_boundaryParticles[j];
				Vector2f r = pi->m_curPosition - pj->m_curPosition;
				if (r.norm() < sl)
					pi->m_neighborList.push_back(pj);
			}
		}
	}
}
void FluidWorld2D::UpdateTimeStepSizeCFL()
{
	float radius = m_fluidParticleRadius;
	float h = m_timeStep;

	// Approximate max. position change due to current velocities
	float maxVel = 0.1f;
	int numParticles = m_numOfParticles;
	float diameter = 2.0f * radius;
	for (int i = 0; i < numParticles; i++)
	{
		Vector2f vel = GetParticle(i)->m_velocity;
		Vector2f accel = GetParticle(i)->m_acceleration;
		float velMag = (vel + accel*h).squaredNorm();
		if (velMag > maxVel)
			maxVel = velMag;
	}

	// boundary particles (if boundary moving)
	/*
	for (unsigned int i = 0; i < m_model->numberOfRigidBodyParticleObjects(); i++)
	{
	FluidModel::RigidBodyParticleObject *rbpo = m_model->getRigidBodyParticleObject(i);
	if (rbpo->m_rigidBody->isDynamic())
	{
	for (unsigned int j = 0; j < rbpo->numberOfParticles(); j++)
	{
	const Vector3r &vel = rbpo->m_v[j];
	const Real velMag = vel.squaredNorm();
	if (velMag > maxVel)
	maxVel = velMag;
	}
	}
	}
	*/

	// Approximate max. time step size 		
	float m_cflFactor = 0.5f;
	float m_cflMaxTimeStepSize = 0.005f;
	float minTimeStepSize = 0.0001f;
	h = m_cflFactor * 0.4f * (diameter / (sqrt(maxVel)));

	h = std::min(h, m_cflMaxTimeStepSize);
	h = std::max(h, minTimeStepSize);

	m_timeStep = h;
}

void FluidWorld2D::StepPBF()
{
	float h = m_timeStep;

	// clear ExternForce
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < m_numOfParticles; i++)
		{
			FParticle2D* pi = m_particles[i];
			pi->m_acceleration = Vector2f(0.0f, -9.8f);
			pi->m_oldPosition = pi->m_curPosition;
			pi->m_oldVelocity = pi->m_velocity;

			pi->m_velocity += pi->m_acceleration * h;
			pi->m_curPosition += pi->m_velocity * h;
			
			pi->m_tempPosition = pi->m_curPosition;
			pi->m_tempVelocity = pi->m_velocity;
		}
	}

	NeighborListUpdate();
	pbfWorld->ConstraintProjection(m_particles, m_boundaryParticles, h, m_kernel);

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < m_numOfParticles; i++)
		{
			m_particles[i]->m_velocity = (m_particles[i]->m_curPosition - m_particles[i]->m_oldPosition) * (1.0f / h);
		}
	}

	//pbfWorld->ComputeXSPHViscosity(m_particles);

	//UpdateTimeStepSizeCFL();

	m_accTimeIntegration += h;
}

void FluidWorld2D::StepPBFonSub1()
{
	float h = m_timeStep;

	// clear ExternForce
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < m_numOfParticles; i++)
		{
			FParticle2D* pi = m_particles[i];
			pi->m_acceleration.setZero();
			pi->m_acceleration[1] = -9.8f;
			pi->m_oldPosition = pi->m_curPosition;
			pi->m_oldVelocity = pi->m_velocity;

			pi->m_velocity += pi->m_acceleration * h;
			pi->m_curPosition += pi->m_velocity * h;
			
			pi->m_tempPosition = pi->m_curPosition;
			pi->m_tempVelocity = pi->m_velocity;
		}
	}

	NeighborListUpdate();
}
void FluidWorld2D::StepPBFonSub2()
{
	float h = m_timeStep;

	pbfWorld->ConstraintProjection(m_particles, m_boundaryParticles, h, m_kernel);

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < m_numOfParticles; i++)
		{
			m_particles[i]->m_velocity = (m_particles[i]->m_curPosition - m_particles[i]->m_oldPosition) * (1.0f / h);
		}
	}

	//pbfWorld->ComputeXSPHViscosity(m_particles);
	//UpdateTimeStepSizeCFL();
	m_accTimeIntegration += h;
}

void FluidWorld2D::StepPBFonSub1WithTF()
{
	float h = m_timeStep;

	// clear ExternForce
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < m_numOfParticles; i++)
		{
			FParticle2D* pi = m_particles[i];
			pi->m_acceleration.setZero();
			pi->m_acceleration[1] = -9.8f;
			pi->m_oldPosition = pi->m_curPosition;
			pi->m_oldVelocity = pi->m_velocity;

			pi->m_velocity += pi->m_acceleration * h;
			pi->m_curPosition += pi->m_velocity * h;

			pi->m_tempPosition = pi->m_curPosition;
			pi->m_tempVelocity = pi->m_velocity;
		}
	}

	NeighborListUpdate();
}
void FluidWorld2D::StepPBFonSub2WithTF()
{
	float h = m_timeStep;

	//pbfWorld->ConstraintProjection(m_particles, m_boundaryParticles, h);

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < m_numOfParticles; i++)
		{
			m_particles[i]->m_velocity = (m_particles[i]->m_curPosition - m_particles[i]->m_oldPosition) * (1.0f / h);
		}
	}

	//pbfWorld->ComputeXSPHViscosity(m_particles);

	//UpdateTimeStepSizeCFL();

	m_accTimeIntegration += h;
}