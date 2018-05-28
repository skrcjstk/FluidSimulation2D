
#include "GL/glew.h"
#include "Visualization\MiniGL.h"
#include "Visualization\Selection.h"
#include "GL/glut.h"
#include "Visualization\PrimitiveBuffer.h"
#include "FluidWorld2D.h"
#include "TimerChrono.h"
#include "PIC2D.h"
#include "OPT2D.h"
#include "PBFC2D.h"
#include <Eigen/Dense>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>

using namespace PBD;
using namespace Eigen;
using namespace std;

void timeStep();
void reset();
void render();
void cleanup();
void buildModel_BreakingDam();
void CreateCoarseBreakingDam(std::vector<Vector2f>& p_damParticles);
void CreateCoarseContainer(std::vector<Vector2f>& p_boundaryParticles);
void CreateFineBreakingDam(std::vector<Vector2f>& p_damParticles);
void CreateFineContainer(std::vector<Vector2f>& p_boundaryParticles);
void AddWall(Vector2f p_min, Vector2f p_max, std::vector<Vector2f>& p_boundaryParticle, float p_particleRadius);

void SavePBFCTrainingData();
void SavePBFCTrainingData2();
void SavePBFCTrainingDataEnv();

TimerChrono timer1;

Primitive spherePrimiCoarse, spherePrimiFine;
GLint context_major_version, context_minor_version;
FluidWorld2D* fineWorld, *coarseWorld;
PIC2D *picForCoarse, *picForFine;
PBFC2D *pbfc2D;
FlowBoundary* wall;
bool doPause = true;

float fineR = 0.025f;
float coarseR = 0.1f;
int fineCoarseRatio = int(coarseR / fineR);
int fineDamWidthCount = 32;
int fineDamHeightCount = 32;
int coarseDamWidthCount = fineDamWidthCount / fineCoarseRatio;
int coarseDamHeightCount = fineDamHeightCount / fineCoarseRatio;
float containerWidth = (coarseDamWidthCount * 4) * coarseR;
float containerHeight = (coarseDamHeightCount * 4) * coarseR;
Vector2f coarseContainerStart, coarseContainerEnd;
Vector2f fineContainerStart, fineContainerEnd;

int frameLimit = 3000;
int accFrameCount = 0;
float accTime = 0.0f;

int maxmaxCountNei = 0;

int main(int argc, char** argv)
{
	//std::srand((unsigned int)std::time(nullptr));
	//printf("OpenMP version: %d\n", _OPENMP);
	
	// OpenGL
	MiniGL::init(argc, argv, 1024, 768, 0, 0, "Fluid demo fine2coarse");
	MiniGL::initLights();
	MiniGL::setClientIdleFunc(50, timeStep);
	MiniGL::setKeyFunc(0, 'r', reset);

	MiniGL::getOpenGLVersion(context_major_version, context_minor_version);

	MiniGL::setClientSceneFunc(render);
	MiniGL::setViewport(40.0f, 0.1f, 500.0f, Vector3f(0.0, 4.5, 13.0), Vector3f(0.0, 4.5, 0.0));

	TwAddVarRW(MiniGL::getTweakBar(), "Pause", TW_TYPE_BOOLCPP, &doPause, " label='Pause' group=Simulation key=SPACE ");

	buildModel_BreakingDam();
	//buildModel_FromFile();
	
	SavePBFCTrainingDataEnv();
	
	if (context_major_version >= 3)
	{
		spherePrimiCoarse.createSphereBuffers((float)coarseR, 8);
		spherePrimiFine.createSphereBuffers((float)fineR, 8);
	}

	glutMainLoop();

	cleanup();

	return 0;
}
void timeStep()
{
	if (doPause)
		return;
	
	//timer1.start();
	//fineWorld->StepPBF();
	//timer1.end();
		
	
	coarseWorld->StepPBF();
	fineWorld->StepPBFonSub1();

	pbfc2D->NeighborBTWTwoResForPBFC(coarseWorld, fineWorld);
	pbfc2D->SolvePBFCConstaints(coarseWorld, fineWorld);
	
	fineWorld->StepPBFonSub2();
	

	//SavePBFCTrainingData();
	SavePBFCTrainingData2();
	
	accFrameCount += 1;
	if (accFrameCount == frameLimit)
	{
		doPause = true;
	}

	accTime += fineWorld->GetTimeStep();
}

void reset()
{
	fineWorld->Reset();
	coarseWorld->Reset();
}
void render()
{
	//MiniGL::coordinateSystem();
	MiniGL::drawTime(accTime);

	float surfaceColor[4] = { 0.2f, 0.2f, 0.2f, 0.1f };
	float kernelColor[4] = { 1.0f, 0.2f, 0.2f, 1.0f };
	float speccolor[4] = { 1.0, 1.0, 1.0, 1.0 };
	float anisotropyColor[4] = { 0.2f, 0.2f, 0.2f, 1.0f };

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, surfaceColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, surfaceColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, speccolor);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);
	glColor3fv(surfaceColor);

	// drawing fine world
	{
		Vector3f translationForFine(-1.75f, 0.0f, 0.0f);
		float fluidColor[4] = { 0.8f, 0.3f, 0.3f, 1.0f };
		for (int i = 0; i < (int)fineWorld->GetNumOfParticles(); i++)
		{
			Vector2f& pos = fineWorld->GetParticle(i)->m_curPosition;
			Vector3f particlePos = Vector3f(pos[0], pos[1], 0.0f) + translationForFine;
			spherePrimiFine.renderSphere(particlePos, fluidColor);
		}

		// drawing fine boundary particles
		float boundaryColor[4] = { 0.2f, 0.2f, 0.2f, 0.5f };
		for (int i = 0; i < fineWorld->GetNumOfBoundaryParticles(); i++)
		{
			Vector2f& pos = fineWorld->GetBoundaryParticle(i)->m_restPosition;
			Vector3f particlePos = Vector3f(pos[0], pos[1], 0.0f) + translationForFine;
			spherePrimiFine.renderSphere(particlePos, boundaryColor);
		}
	}

	// drawing coarse world
	{
		Vector3f translationForCoarse(1.75f, 0.0f, 0.0f);
		float fluidColorForCoarse[4] = { 0.0f, 0.2f, 0.7f, 1.0f };
		for (int i = 0; i < coarseWorld->GetNumOfParticles(); i++)
		{
			Vector2f& pos = coarseWorld->GetParticle(i)->m_curPosition;
			Vector3f particlePos = Vector3f(pos[0], pos[1], 0.0f) + translationForCoarse;
			spherePrimiCoarse.renderSphere(particlePos, fluidColorForCoarse);
		}

		// drawing coarse boundary particles
		float boundaryColor[4] = { 0.2f, 0.2f, 0.2f, 0.5f };
		for (int i = 0; i < coarseWorld->GetNumOfBoundaryParticles(); i++)
		{
			Vector2f& pos = coarseWorld->GetBoundaryParticle(i)->m_restPosition;
			Vector3f particlePos = Vector3f(pos[0], pos[1], 0.0f) + translationForCoarse;
			spherePrimiCoarse.renderSphere(particlePos, boundaryColor);
		}
	}

	// drawing union world
	{
		Vector3f translationForUnion(0.0f, 3.5f, 0.0f);

		float fluidColor[4] = { 0.8f, 0.3f, 0.3f, 1.0f };
		for (int i = 0; i < (int)fineWorld->GetNumOfParticles(); i++)
		{
			Vector2f& pos = fineWorld->GetParticle(i)->m_curPosition;
			Vector3f particlePos = Vector3f(pos[0], pos[1], 0.0f) + translationForUnion;
			spherePrimiFine.renderSphere(particlePos, fluidColor);
		}

		float fluidColorForCoarse[4] = { 0.0f, 0.2f, 0.7f, 1.0f };
		for (int i = 0; i < coarseWorld->GetNumOfParticles(); i++)
		{
			Vector2f& pos = coarseWorld->GetParticle(i)->m_curPosition;
			Vector3f particlePos = Vector3f(pos[0], pos[1], 0.0f) + translationForUnion;
			spherePrimiCoarse.renderSphere(particlePos, fluidColorForCoarse);
		}

		// drawing fine boundary particles
		float boundaryColor[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
		for (int i = 0; i < fineWorld->GetNumOfBoundaryParticles(); i++)
		{
			Vector2f& pos = fineWorld->GetBoundaryParticle(i)->m_restPosition;
			Vector3f particlePos = Vector3f(pos[0], pos[1], 0.0f) + translationForUnion;
			spherePrimiFine.renderSphere(particlePos, boundaryColor);
		}

		// drawing fine boundary particles
		float boundaryColorForCoarse[4] = { 0.2f, 0.2f, 0.2f, 0.5f };
		for (int i = 0; i < coarseWorld->GetNumOfBoundaryParticles(); i++)
		{
			Vector2f& pos = coarseWorld->GetBoundaryParticle(i)->m_restPosition;
			Vector3f particlePos = Vector3f(pos[0], pos[1], 0.0f) + translationForUnion;
			spherePrimiCoarse.renderSphere(particlePos, boundaryColorForCoarse);
		}
	}

}
void cleanup()
{
	if (context_major_version >= 3)
	{
		spherePrimiCoarse.releaseBuffers();
		spherePrimiFine.releaseBuffers();
	}
	timer1.printAvg("fineWorld->StepPBF");
}

void buildModel_BreakingDam()
{
	float timeStep = 0.0025f;

	// sub domain creation
	std::vector<Vector2f> fineBoundaryParticles;
	std::vector<Vector2f> fineFluidParticles;
	CreateFineBreakingDam(fineFluidParticles);
	CreateFineContainer(fineBoundaryParticles);

	fineWorld = new FluidWorld2D();
	fineWorld->SetTimeStep(timeStep);
	fineWorld->CreateBoundaryParticles(fineBoundaryParticles, coarseR);
	fineWorld->CreateFluidParticles(fineFluidParticles, fineR);

	// main domain creation
	std::vector<Vector2f> coarseBoundaryParticles;
	std::vector<Vector2f> coarseFluidParticles;
	CreateCoarseBreakingDam(coarseFluidParticles);
	CreateCoarseContainer(coarseBoundaryParticles);
	
	coarseWorld = new FluidWorld2D();
	coarseWorld->SetTimeStep(timeStep);
	coarseWorld->CreateBoundaryParticles(coarseBoundaryParticles, coarseR);
	coarseWorld->CreateFluidParticles(coarseFluidParticles, coarseR);

	pbfc2D = new PBFC2D(coarseWorld, fineWorld);
	pbfc2D->SetSearchRange(2.0f * coarseWorld->GetParticleRadius());
}

void CreateCoarseBreakingDam(std::vector<Vector2f>& p_damParticles)
{
	p_damParticles.resize(coarseDamWidthCount*coarseDamHeightCount);

	float diam = 2.0f * coarseR;
	float startX = -0.5f * containerWidth + 0.2f * containerHeight;
	startX += coarseR;
	float startY = diam + diam + diam;
	startY += coarseR;

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int j = 0; j < coarseDamHeightCount; j++)
		{
			for (int i = 0; i < coarseDamWidthCount; i++)
			{
				p_damParticles[j*coarseDamWidthCount + i] = diam*Eigen::Vector2f((float)i, (float)j) + Eigen::Vector2f(startX, startY);
			}
		}
	}
}
void CreateCoarseContainer(std::vector<Vector2f>& p_boundaryParticles)
{
	float x1 = -containerWidth / 2.0f;
	float x2 = containerWidth / 2.0f;
	float y1 = 0.0f;
	float y2 = containerHeight;

	coarseContainerStart[0] = x1;
	coarseContainerStart[1] = y1;
	coarseContainerEnd[0] = x2;
	coarseContainerEnd[1] = y2;

	// Floor
	AddWall(Vector2f(x1, y1), Vector2f(x2, y1), p_boundaryParticles, coarseR);
	// Top
	//AddWall(Vector2f(x1, y2), Vector2f(x2, y2), p_boundaryParticles, coarseR);
	// Left
	AddWall(Vector2f(x1, y1), Vector2f(x1, y2), p_boundaryParticles, coarseR);
	// Right
	AddWall(Vector2f(x2, y1), Vector2f(x2, y2), p_boundaryParticles, coarseR);
}
void CreateFineBreakingDam(std::vector<Vector2f>& p_damParticles)
{
	p_damParticles.resize(fineDamWidthCount*fineDamHeightCount);

	float diam = 2.0f * fineR;
	float coarseDiam = 2.0f * coarseR;
	float startX = -0.5f * containerWidth + 0.2f * containerHeight;
	startX += fineR;
	float startY = coarseDiam + coarseDiam + coarseDiam;
	startY += fineR;

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int j = 0; j < fineDamHeightCount; j++)
		{
			for (int i = 0; i < fineDamWidthCount; i++)
			{
				p_damParticles[j*fineDamWidthCount + i] = diam*Eigen::Vector2f((float)i, (float)j) + Eigen::Vector2f(startX, startY);
			}
		}
	}
}
void CreateFineContainer(std::vector<Vector2f>& p_boundaryParticles)
{
	float x1 = -containerWidth / 2.0f;
	float x2 = containerWidth / 2.0f;
	x1 += fineCoarseRatio * 0.5f * fineR;
	x2 -= fineCoarseRatio * 0.5f * fineR;
	float y1 = 0.0f;
	float y2 = containerHeight;
	y1 += fineCoarseRatio * 0.5f * fineR;
	y2 -= fineCoarseRatio * 0.5f * fineR;

	fineContainerStart[0] = x1;
	fineContainerStart[1] = y1;
	fineContainerEnd[0] = x2;
	fineContainerEnd[1] = y2;

	// Floor
	AddWall(Vector2f(x1, y1), Vector2f(x2, y1), p_boundaryParticles, fineR);
	// Top
	//AddWall(Vector3f(x1, y2, z1), Vector3f(x2, y2, z2), p_boundaryParticles, fineR);
	// Left
	AddWall(Vector2f(x1, y1), Vector2f(x1, y2), p_boundaryParticles, fineR);
	// Right
	AddWall(Vector2f(x2, y1), Vector2f(x2, y2), p_boundaryParticles, fineR);
}
void AddWall(Vector2f p_min, Vector2f p_max, std::vector<Vector2f>& p_boundaryParticle, float p_particleRadius)
{
	Vector2f diff = p_max - p_min;
	float diameter = 2 * p_particleRadius;

	int countX = (int)(diff[0] / diameter) + 1;
	int countY = (int)(diff[1] / diameter) + 1;

	int startIndex = (int)p_boundaryParticle.size();
	p_boundaryParticle.resize(startIndex + countX*countY);

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < countX; i++)
		{
			for (int j = 0; j < countY; j++)
			{
				const Vector2f position = p_min + Vector2f(i*diameter, j*diameter);
				p_boundaryParticle[startIndex + i*countY + j] = position;
			}
		}
	}
}

void SavePBFCTrainingDataEnv()
{
	float fbuf[2];
	int ibuf[1];

	FILE* fpEnvF = fopen("./PBFC2D_SD15/BoundaryEnv.dat", "wb");

	ibuf[0] = fineWorld->GetNumOfBoundaryParticles();
	fwrite(ibuf, sizeof(int), 1, fpEnvF);
	for (int i = 0; i < fineWorld->GetNumOfBoundaryParticles(); i++)
	{
		fbuf[0] = fineWorld->GetBoundaryParticle(i)->m_mass;
		fwrite(fbuf, sizeof(float), 1, fpEnvF);
		fbuf[0] = fineWorld->GetBoundaryParticle(i)->m_curPosition[0];
		fbuf[1] = fineWorld->GetBoundaryParticle(i)->m_curPosition[1];
		fwrite(fbuf, sizeof(float), 2, fpEnvF);
	}
	fclose(fpEnvF);
	printf("EnvForFine are saved.\n");
}

void SavePBFCTrainingData()
{
	std::vector<FParticle2D*>& fineP = fineWorld->GetParticleList();
	std::vector<FParticle2D*>& coarseP = coarseWorld->GetParticleList();

	string frameIdx = std::to_string(accFrameCount) + ".dat";
	string NeighborIdxPath = "./PBFC2D_SD15/NeighborInfo_";
	string FineParticleInfoPath = "./PBFC2D_SD15/FineParticleInfo_";
	string CoarseParticleInfoPath = "./PBFC2D_SD15/CoarseParticleInfo_";
	string GTPosPath = "./PBFC2D_SD15/GTPos_";
	string GTVelPath = "./PBFC2D_SD15/GTVel_";

	FILE* fpNeiIdx = fopen((NeighborIdxPath + frameIdx).c_str(), "wb");
	FILE* fpFineP = fopen((FineParticleInfoPath + frameIdx).c_str(), "wb");
	FILE* fpCoarseP = fopen((CoarseParticleInfoPath + frameIdx).c_str(), "wb");
	FILE* fpGTPos = fopen((GTPosPath + frameIdx).c_str(), "wb");
	FILE* fpGTVel = fopen((GTVelPath + frameIdx).c_str(), "wb");

	float fbuf[2];
	int ibuf[1];

	// num of fine particles 
	ibuf[0] = (int)fineP.size();
	fwrite(ibuf, sizeof(int), 1, fpFineP);
	fwrite(ibuf, sizeof(int), 1, fpGTPos);
	fwrite(ibuf, sizeof(int), 1, fpGTVel);

	int nNei, maxCountNei = 0;
	int nNeiFP, nNeiBP;
	// save fineP's information
	for (int i = 0; i < (int)fineP.size(); i++)
	{
		FParticle2D* pi = fineP[i];

		// particle information save (mass, tempPosition, tempVelocity)
		fbuf[0] = pi->m_mass;
		fwrite(fbuf, sizeof(float), 1, fpFineP);
		fbuf[0] = pi->m_tempPosition[0];	fbuf[1] = pi->m_tempPosition[1];
		fwrite(fbuf, sizeof(float), 2, fpFineP);
		fbuf[0] = pi->m_tempVelocity[0];	fbuf[1] = pi->m_tempVelocity[1];
		fwrite(fbuf, sizeof(float), 2, fpFineP);

		// neighbor particle info (neighborList Size, neighbgor(coarse) pidx)
		nNeiFP = (int)pbfc2D->m_neighListwithCoarseFP[i].size();
		nNeiBP = (int)pbfc2D->m_neighListwithCoarseBP[i].size();
		nNei = ibuf[0] = nNeiFP + nNeiBP;
		fwrite(ibuf, sizeof(int), 1, fpNeiIdx);
		for (int j = 0; j < nNeiFP; j++)
		{
			ibuf[0] = (int)Pid::Fluid; // Fluid = 0
			ibuf[1] = pbfc2D->m_neighListwithCoarseFP[i][j];
			fwrite(ibuf, sizeof(int), 2, fpNeiIdx);
		}
		for (int j = 0; j < nNeiBP; j++)
		{
			ibuf[0] = (int)Pid::Boundary; // Boundary = 1
			ibuf[1] = pbfc2D->m_neighListwithCoarseBP[i][j];
			fwrite(ibuf, sizeof(int), 2, fpNeiIdx);
		}

		maxCountNei < nNei ? maxCountNei = nNei : maxCountNei = maxCountNei;

		// GTPos & GTVel
		fbuf[0] = pi->m_curPosition[0] - pi->m_tempPosition[0];
		fbuf[1] = pi->m_curPosition[1] - pi->m_tempPosition[1];
		fwrite(fbuf, sizeof(float), 2, fpGTPos);
		fbuf[0] = fbuf[0] * (1.0f / fineWorld->GetTimeStep());
		fbuf[1] = fbuf[1] * (1.0f / fineWorld->GetTimeStep());
		fwrite(fbuf, sizeof(float), 2, fpGTVel);
	}

	ibuf[0] = (int)coarseP.size();
	fwrite(ibuf, sizeof(int), 1, fpCoarseP);
	// save coarseP's information
	for (int i = 0; i < (int)coarseP.size(); i++)
	{
		FParticle2D* pi = coarseP[i];
		
		// particle information save (mass, tempPosition, tempVelocity)
		fbuf[0] = pi->m_mass;
		fwrite(fbuf, sizeof(float), 1, fpCoarseP);
		fbuf[0] = pi->m_curPosition[0];	fbuf[1] = pi->m_curPosition[1];
		fwrite(fbuf, sizeof(float), 2, fpCoarseP);
		fbuf[0] = pi->m_velocity[0];	fbuf[1] = pi->m_velocity[1];
		fwrite(fbuf, sizeof(float), 2, fpCoarseP);
	}
	
	maxmaxCountNei < maxCountNei ? maxmaxCountNei = maxCountNei : maxmaxCountNei = maxmaxCountNei;
	printf("(%d)frame maxmaxCountNei: %d \n", accFrameCount, maxmaxCountNei);

	fclose(fpNeiIdx);
	fclose(fpFineP);
	fclose(fpCoarseP);
	fclose(fpGTPos);
	fclose(fpGTVel);
}

void SavePBFCTrainingData2()
{
	std::vector<FParticle2D*>& fineP = fineWorld->GetParticleList();
	std::vector<FParticle2D*>& coarseP = coarseWorld->GetParticleList();

	string frameIdx = std::to_string(accFrameCount) + ".dat";
	string FineParticleInfoPath = "./PBFC2D_SD15/FineParticleInfo_";
	string CoarseParticleInfoPath = "./PBFC2D_SD15/CoarseParticleInfo_";

	FILE* fpFineP = fopen((FineParticleInfoPath + frameIdx).c_str(), "wb");
	FILE* fpCoarseP = fopen((CoarseParticleInfoPath + frameIdx).c_str(), "wb");
	
	float fbuf[2];
	int ibuf[1];

	// num of fine particles 
	ibuf[0] = (int)fineP.size();
	fwrite(ibuf, sizeof(int), 1, fpFineP);

	// save fineP's information
	for (int i = 0; i < (int)fineP.size(); i++)
	{
		FParticle2D* pi = fineP[i];

		// particle information save (mass, tempPosition, tempVelocity)
		fbuf[0] = pi->m_mass;
		fwrite(fbuf, sizeof(float), 1, fpFineP);
		fbuf[0] = pi->m_curPosition[0];	fbuf[1] = pi->m_curPosition[1];
		fwrite(fbuf, sizeof(float), 2, fpFineP);
		fbuf[0] = pi->m_tempVelocity[0];	fbuf[1] = pi->m_tempVelocity[1];
		fwrite(fbuf, sizeof(float), 2, fpFineP);
	}

	ibuf[0] = (int)coarseP.size();
	fwrite(ibuf, sizeof(int), 1, fpCoarseP);
	// save coarseP's information
	for (int i = 0; i < (int)coarseP.size(); i++)
	{
		FParticle2D* pi = coarseP[i];

		// particle information save (mass, tempPosition, tempVelocity)
		fbuf[0] = pi->m_mass;
		fwrite(fbuf, sizeof(float), 1, fpCoarseP);
		fbuf[0] = pi->m_curPosition[0];	fbuf[1] = pi->m_curPosition[1];
		fwrite(fbuf, sizeof(float), 2, fpCoarseP);
		fbuf[0] = pi->m_velocity[0];	fbuf[1] = pi->m_velocity[1];
		fwrite(fbuf, sizeof(float), 2, fpCoarseP);
	}

	fclose(fpFineP);
	fclose(fpCoarseP);

	printf("%d frame saved.\n", accFrameCount);
}
