/*
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
void buildModel_FromFile();
void CreateCoarseBreakingDam(std::vector<Vector2f>& p_damParticles);
void CreateCoarseContainer(std::vector<Vector2f>& p_boundaryParticles);
void CreateFineBreakingDam(std::vector<Vector2f>& p_damParticles);
void CreateFineContainer(std::vector<Vector2f>& p_boundaryParticles);
void AddWall(Vector2f p_min, Vector2f p_max, std::vector<Vector2f>& p_boundaryParticle, float p_particleRadius);

void LoadFineFluidInfo(int idx, std::vector<Vector2f>& p_fineBoundary, 
	std::vector<Vector2f>& p_fineFluid, std::vector<Vector2f>& p_fineFluidVel, std::vector<Vector2f>& p_coarseBoundary);
void SaveFluidEnvInfo();
void SaveFineFluidInfo(int idx);

void ParticleTrainingDataSave();
void BoundaryParticleInfoSave();
void ParticleTrainingDataSave_fineOnly_simulBefore();
void ParticleTrainingDataSave_fineOnly_simulAfter();

Primitive spherePrimiCoarse, spherePrimiFine;
GLint context_major_version, context_minor_version;
FluidWorld2D* fineWorld, *coarseWorld;
PIC2D *picForCoarse, *picForFine;
PBFC2D *pbfc2D;
FlowBoundary* wall;
bool doPause = true;

float fineR = 0.025f;
float coarseR = 0.05f;
int fineDamWidth = 30;
int fineDamHeight = 30;
int coarseDamWidth = fineDamWidth / 2;
int coarseDamHeight = fineDamHeight / 2;
float containerWidth = (coarseDamWidth * 4) * coarseR;
float containerHeight = (coarseDamWidth * 4) * coarseR;
Vector2f coarseContainerStart, coarseContainerEnd;
Vector2f fineContainerStart, fineContainerEnd;

float gDxForF;
Vector2f gStartForF, gEndForF, gSizeForF;
float gDxForC;
Vector2f gStartForC, gEndForC, gSizeForC;

int frameLimit = 500;
int accFrameCount = 0;

int maxmaxF=0, maxmaxB=0, maxmaxFB=0;


int main(int argc, char** argv)
{
	std::srand((unsigned int)std::time(nullptr));
	printf("OpenMP version: %d\n", _OPENMP);

	// OpenGL 
	MiniGL::init(argc, argv, 1024, 768, 0, 0, "Fluid demo");
	MiniGL::initLights();
	MiniGL::setClientIdleFunc(50, timeStep);
	MiniGL::setKeyFunc(0, 'r', reset);

	MiniGL::getOpenGLVersion(context_major_version, context_minor_version);

	MiniGL::setClientSceneFunc(render);
	MiniGL::setViewport(40.0f, 0.1f, 500.0f, Vector3f(0.0, 4.0, 15.0), Vector3f(0.0, 4.0, 0.0));

	TwAddVarRW(MiniGL::getTweakBar(), "Pause", TW_TYPE_BOOLCPP, &doPause, " label='Pause' group=Simulation key=SPACE ");

	buildModel_BreakingDam();
	//buildModel_FromFile();

	if (context_major_version >= 3)
	{
		spherePrimiCoarse.createSphereBuffers((float)coarseR, 8);
		spherePrimiFine.createSphereBuffers((float)fineR, 8);
	}

	BoundaryParticleInfoSave();

	glutMainLoop();

	cleanup();
		
	return 0;
}


void timeStep()
{
	if (doPause)
		return;

	for (unsigned int i = 0; i < 1; i++)
	{
		fineWorld->StepPBF();
		//picForFine->AssignCells(fineWorld);
		//picForFine->Map_P2G(fineWorld);

		//coarseWorld->StepPBFonSub1();
		//pbfc2D->NeighborBTWTwoResForPBFC(fineWorld, coarseWorld);
		//pbfc2D->SolvePBFCConstaints(fineWorld, coarseWorld);
		//coarseWorld->StepPBFonSub2();
		
		//picForCoarse->AssignCells(coarseWorld);
		//picForCoarse->Map_P2G(coarseWorld);
	}
	//doPause = !doPause;

	// training data (particle info) save
	//pbfc2D->UpdateTrainingData(fineWorld, coarseWorld);
	//ParticleTrainingDataSave();

	//ParticleTrainingDataSave_fineOnly_simulBefore();
	//ParticleTrainingDataSave_fineOnly_simulAfter();

	accFrameCount += 1;
	if (accFrameCount == frameLimit)
	{
		//SaveFineFluidInfo(accFrameCount);
		doPause = true;
	}
}

void reset()
{
	fineWorld->Reset();
	coarseWorld->Reset();
}
void render()
{
	//MiniGL::coordinateSystem();
	MiniGL::drawTime(0.0f);

	float surfaceColor[4] = { 0.2f, 0.2f, 0.2f, 0.1f };
	float kernelColor[4] = { 1.0f, 0.2f, 0.2f, 1.0f };
	float speccolor[4] = { 1.0, 1.0, 1.0, 1.0 };
	float anisotropyColor[4] = { 0.2f, 0.2f, 0.2f, 1.0f };

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, surfaceColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, surfaceColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, speccolor);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);
	glColor3fv(surfaceColor);

	// drawing fine fluid world
	{
		glPointSize(6.0f);
		Vector2f translationForFine(-0.0f, 0.0f);

		glPushMatrix();
		glTranslatef(translationForFine[0], translationForFine[1], 0.0f);
		glBegin(GL_POINTS);
		float fluidColor[4] = { 0.0f, 0.7f, 0.7f, 0.2f };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, fluidColor);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, fluidColor);
		for (int i = 0; i < fineWorld->GetNumOfParticles(); i++)
		{
			Vector2f& pos = fineWorld->GetParticle(i)->m_curPosition;
			glVertex2f(pos[0], pos[1]);
		}
		glEnd();
		glPopMatrix();

		// drawing fine boundary particles
		glPushMatrix();
		glTranslatef(translationForFine[0], translationForFine[1], 0.0f);
		glBegin(GL_POINTS);
		float boundaryColor[4] = { 0.2f, 0.2f, 0.2f, 0.5f };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, boundaryColor);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, boundaryColor);
		for (int i = 0; i < fineWorld->GetNumOfBoundaryParticles(); i++)
		{
			Vector2f& pos = fineWorld->GetBoundaryParticle(i)->m_restPosition;
			glVertex2f(pos[0], pos[1]);
		}
		glEnd();
		glPopMatrix();


		// drawing grid cell
		glPushMatrix();
		glTranslatef(translationForFine[0], translationForFine[1], 0.0f);
		glBegin(GL_LINES);
		float lineColor[4] = { 0.2f, 0.2f, 0.8f, 0.8f };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, lineColor);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, lineColor);
		int xCount = (int)(gSizeForF[0] / gDxForF) + 1;
		int yCount = (int)(gSizeForF[1] / gDxForF) + 1;
		for (int i = 0; i < xCount; i++)
		{
			Vector2f from(gStartForF[0] + gDxForF * i, gStartForF[1]);
			Vector2f end(from[0], gEndForF[1]);
			glVertex2f(from[0], from[1]);
			glVertex2f(end[0], end[1]);
		}
		for (int i = 0; i < yCount; i++)
		{
			Vector2f from(gStartForF[0], gStartForF[1] + gDxForF * i);
			Vector2f end(gEndForF[0], from[1]);
			glVertex2f(from[0], from[1]);
			glVertex2f(end[0], end[1]);
		}
		glEnd();
		glPopMatrix();

		// render arrow
		Vector2f dxy = picForFine->GetDxDy();
		float head_len = 0.3f*dxy[0];
		int ng = picForFine->GetGridSize();
		for (int g = 0; g < ng; g++)
		{
			Vector2f pos = picForFine->GetGridPos(g) + translationForFine;
			Vector2f end = (pos + 0.05f * picForFine->GetVelocity(g));
			spherePrimiCoarse.renderArrow2D(pos , end, head_len);
		}
	}
	
	// drawing coarse fluid world
	{
		glPointSize(10.0f);
		Vector2f translationForCoarse(0.0f, 0.0f);
		
		glPushMatrix();
		glTranslatef(translationForCoarse[0], translationForCoarse[1], 0.0f);
		glBegin(GL_POINTS);
		float fluidColorForCoarse[4] = { 0.0f, 0.2f, 0.7f, 0.5f };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, fluidColorForCoarse);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, fluidColorForCoarse);
		for (int i = 0; i < coarseWorld->GetNumOfParticles(); i++)
		{
			Vector2f& pos = coarseWorld->GetParticle(i)->m_curPosition;
			glVertex2f(pos[0], pos[1]);
		}
		glEnd();
		glPopMatrix();

		// drawing coarse boundary particles
		glPushMatrix();
		glTranslatef(translationForCoarse[0], translationForCoarse[1], 0.0f);
		glBegin(GL_POINTS);
		float boundaryColorForCoarse[4] = { 0.2f, 0.2f, 0.2f, 0.5f };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, boundaryColorForCoarse);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, boundaryColorForCoarse);
		for (int i = 0; i < coarseWorld->GetNumOfBoundaryParticles(); i++)
		{
			Vector2f& pos = coarseWorld->GetBoundaryParticle(i)->m_curPosition;
			glVertex2f(pos[0], pos[1]);
		}
		glEnd();
		glPopMatrix();

		// drawing OPT2D & PIC grid cell
		glPushMatrix();
		glTranslatef(translationForCoarse[0], translationForCoarse[1], 0.0f);
		glBegin(GL_LINES);
		float lineColor[4] = { 0.2f, 0.2f, 0.8f, 0.8f };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, lineColor);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, lineColor);
		int xCount = (int)(gSizeForC[0] / gDxForC) + 1;
		int yCount = (int)(gSizeForC[1] / gDxForC) + 1;
		for (int i = 0; i < xCount; i++)
		{
			Vector2f from(gStartForC[0] + gDxForC * i, gStartForC[1]);
			Vector2f end(from[0], gEndForC[1]);
			glVertex2f(from[0], from[1]);
			glVertex2f(end[0], end[1]);
		}
		for (int i = 0; i < yCount; i++)
		{
			Vector2f from(gStartForC[0], gStartForC[1] + gDxForC * i);
			Vector2f end(gEndForC[0], from[1]);
			glVertex2f(from[0], from[1]);
			glVertex2f(end[0], end[1]);
		}
		glEnd();
		glPopMatrix();

		// render arrow
		Vector2f dxy = picForCoarse->GetDxDy();
		float head_len = 0.3f*dxy[0];
		int ng = picForCoarse->GetGridSize();
		for (int g = 0; g < ng; g++)
		{
			Vector2f pos = picForCoarse->GetGridPos(g) + translationForCoarse;
			Vector2f end = (pos + 0.05f * picForCoarse->GetVelocity(g));
			spherePrimiCoarse.renderArrow2D(pos , end, head_len);
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
}
void buildModel_BreakingDam()
{
	float timeStep = 0.0025f;

	// main domain creation
	std::vector<Vector2f> boundaryParticles;
	std::vector<Vector2f> damParticles;
	CreateFineBreakingDam(damParticles);
	CreateFineContainer(boundaryParticles);

	fineWorld = new FluidWorld2D();
	fineWorld->SetTimeStep(timeStep);
	fineWorld->CreateBoundaryParticles(boundaryParticles, fineR);
	fineWorld->CreateFluidParticles(damParticles, fineR);
	
	// PIC grid for fine grid
	gDxForF = 2.0f * fineR;
	gStartForF = fineContainerStart - Vector2f(6.0f * gDxForF, 6.0f * gDxForF);
	gEndForF = fineContainerEnd + Vector2f(6.0f*gDxForF, 6.0f*gDxForF);
	gSizeForF = gEndForF - gStartForF;

	picForFine = new PIC2D();
	picForFine->Initialize(fineWorld, gStartForF, gSizeForF, Vector2i((int)(gSizeForF[0] / gDxForF), (int)(gSizeForF[1] / gDxForF)), 1.0);
	picForFine->AssignBoundary(fineWorld->GetBoundaryParticleList());

	// sub-domain creation
	std::vector<Vector2f> subBoundary;
	std::vector<Vector2f> a;
	CreateCoarseContainer(subBoundary);
	
	// create sub-domain fluid particles using OPT2D
	gDxForC = 2.0f * coarseR;
	gStartForC = coarseContainerStart - Vector2f(3.0f * gDxForC, 3.0f * gDxForC);
	gEndForC = coarseContainerEnd + Vector2f(3.0f*gDxForC, 3.0f*gDxForC);
	gSizeForC = gEndForC - gStartForC;

	wall = new FlowBoundary();
	wall->SetFluidThreshold(0.5f);
	wall->CreateFlowBoundary(gStartForC, gEndForC, gDxForC);
	wall->CheckInBoundary(subBoundary);
	wall->CreateCoarsePS(fineWorld, a);

	coarseWorld = new FluidWorld2D();
	coarseWorld->SetTimeStep(timeStep);
	coarseWorld->CreateFluidParticles(a, coarseR);
	coarseWorld->CreateBoundaryParticles(subBoundary, coarseR);

	wall->NeighborSearchBTWTwoRes(fineWorld, coarseWorld);
	wall->InterpolateVelocity(fineWorld, coarseWorld);

	picForCoarse = new PIC2D();
	picForCoarse->Initialize(coarseWorld, gStartForC, gSizeForC, Vector2i((int)(gSizeForC[0] / gDxForC), (int)(gSizeForC[1] / gDxForC)), 1.0);
	picForCoarse->AssignBoundary(coarseWorld->GetBoundaryParticleList());
	
	// Save FluidEnv(Boundary particles) 
	//SaveFluidEnvInfo();

	pbfc2D = new PBFC2D(fineWorld, coarseWorld);

}

void buildModel_FromFile()
{
	float timeStep = 0.0025f;

	// main domain creation
	std::vector<Vector2f> boundaryParticles;
	std::vector<Vector2f> damParticles;
	std::vector<Vector2f> damParticlesVel;

	std::vector<Vector2f> subBoundary;
	std::vector<Vector2f> a;

	LoadFineFluidInfo(120, boundaryParticles, damParticles, damParticlesVel, subBoundary);

	fineWorld = new FluidWorld2D();
	fineWorld->SetTimeStep(timeStep);
	fineWorld->CreateBoundaryParticles(boundaryParticles, fineR);
	fineWorld->CreateFluidParticles(damParticles, fineR);
	
	// velocity 
	for (int i = 0; i < fineWorld->GetNumOfParticles(); i++)
		fineWorld->GetParticle(i)->m_velocity = damParticlesVel[i];

	// PIC grid for coarse grid
	gDxForF = 2.0f * fineR;
	gStartForF = fineContainerStart - Vector2f(6.0f * gDxForF, 6.0f * gDxForF);
	gEndForF = fineContainerEnd + Vector2f(6.0f*gDxForF, 6.0f*gDxForF);
	gSizeForF = gEndForF - gStartForF;

	picForFine = new PIC2D();
	picForFine->Initialize(fineWorld, gStartForF, gSizeForF, Vector2i((int)(gSizeForF[0] / gDxForF), (int)(gSizeForF[1] / gDxForF)), 1.0);
	picForFine->AssignBoundary(fineWorld->GetBoundaryParticleList());

	// create sub-domain fluid particles using OPT2D
	gDxForC = 2.0f * coarseR;
	gStartForC = coarseContainerStart - Vector2f(3.0f * gDxForC, 3.0f * gDxForC);
	gEndForC = coarseContainerEnd + Vector2f(3.0f*gDxForC, 3.0f*gDxForC);
	gSizeForC = gEndForC - gStartForC;

	wall = new FlowBoundary();
	wall->SetFluidThreshold(0.5f);
	wall->CreateFlowBoundary(gStartForC, gEndForC, gDxForC);
	wall->CheckInBoundary(subBoundary);
	wall->CreateCoarsePS(fineWorld, a);

	coarseWorld = new FluidWorld2D();
	coarseWorld->SetTimeStep(timeStep);
	coarseWorld->CreateFluidParticles(a, coarseR);
	coarseWorld->CreateBoundaryParticles(subBoundary, coarseR);

	wall->NeighborSearchBTWTwoRes(fineWorld, coarseWorld);
	wall->InterpolateVelocity(fineWorld, coarseWorld);

	picForCoarse = new PIC2D();
	picForCoarse->Initialize(coarseWorld, gStartForC, gSizeForC, Vector2i((int)(gSizeForC[0] / gDxForC), (int)(gSizeForC[1] / gDxForC)), 1.0);
	picForCoarse->AssignBoundary(coarseWorld->GetBoundaryParticleList());
}

void LoadFineFluidInfo(int idx, std::vector<Vector2f>& p_fineBoundary, std::vector<Vector2f>& p_fineFluid, std::vector<Vector2f>& p_fineFluidVel, std::vector<Vector2f>& p_coarseBoundary)
{
	FILE* fpFluidF, *fpEnvC, *fpEnvF;
	string filename = "FluidForFine_" + std::to_string(idx) + ".dat";
	float fbuf[2];
	int ibuf[1];

	fpFluidF = fopen(filename.c_str(), "rb");
	fpEnvF = fopen("EnvForFine.dat", "rb");
	fpEnvC = fopen("EnvForCoarse.dat", "rb");

	fread(ibuf, sizeof(int), 1, fpEnvF);
	fread(fbuf, sizeof(float), 2, fpEnvF);
	fineContainerStart[0] = fbuf[0];	fineContainerStart[1] = fbuf[1];
	fread(fbuf, sizeof(float), 2, fpEnvF);
	fineContainerEnd[0] = fbuf[0];	fineContainerEnd[1] = fbuf[1];
	for (int i = 0; i < ibuf[0]; i++)
	{
		fread(fbuf, sizeof(float), 2, fpEnvF);
		p_fineBoundary.push_back(Vector2f(fbuf[0], fbuf[1]));
	}

	fread(ibuf, sizeof(int), 1, fpFluidF);
	for (int i = 0; i < ibuf[0]; i++)
	{
		fread(fbuf, sizeof(float), 2, fpFluidF);
		p_fineFluid.push_back(Vector2f(fbuf[0], fbuf[1]));
		fread(fbuf, sizeof(float), 2, fpFluidF);
		p_fineFluidVel.push_back(Vector2f(fbuf[0], fbuf[1]));
	}

	fread(ibuf, sizeof(int), 1, fpEnvC);
	fread(fbuf, sizeof(float), 2, fpEnvC);
	coarseContainerStart[0] = fbuf[0];	coarseContainerStart[1] = fbuf[1];
	fread(fbuf, sizeof(float), 2, fpEnvC);
	coarseContainerEnd[0] = fbuf[0];	coarseContainerEnd[1] = fbuf[1];
	for (int i = 0; i < ibuf[0]; i++)
	{
		fread(fbuf, sizeof(float), 2, fpEnvC);
		p_coarseBoundary.push_back(Vector2f(fbuf[0], fbuf[1]));
	}

	fclose(fpFluidF);
	fclose(fpEnvF);
	fclose(fpEnvC);

}
void SaveFineFluidInfo(int idx)
{
	float fbuf[2];
	int ibuf[1];
	FILE* fpFluidF;
	string filename = "FluidForFine_" + std::to_string(idx) + ".dat";
	
	fpFluidF = fopen(filename.c_str(), "wb");
	ibuf[0] = fineWorld->GetNumOfParticles();
	fwrite(ibuf, sizeof(int), 1, fpFluidF);
	for (int i = 0; i < ibuf[0]; i++)
	{
		fbuf[0] = fineWorld->GetParticle(i)->m_curPosition[0];
		fbuf[1] = fineWorld->GetParticle(i)->m_curPosition[1];
		fwrite(fbuf, sizeof(float), 2, fpFluidF);
		fbuf[0] = fineWorld->GetParticle(i)->m_velocity[0];
		fbuf[1] = fineWorld->GetParticle(i)->m_velocity[1];
		fwrite(fbuf, sizeof(float), 2, fpFluidF);
	}
	fclose(fpFluidF);
	printf("FluidForFine_%d is saved.\n", idx);
}
void SaveFluidEnvInfo()
{
	float fbuf[2];
	int ibuf[1];
	FILE* fpEnvC, *fpEnvF;

	fpEnvC = fopen("EnvForCoarse.dat", "wb");
	fpEnvF = fopen("EnvForFine.dat", "wb");

	ibuf[0] = fineWorld->GetNumOfBoundaryParticles();
	fwrite(ibuf, sizeof(int), 1, fpEnvF);
	fbuf[0] = fineContainerStart[0];
	fbuf[1] = fineContainerStart[1];
	fwrite(fbuf, sizeof(float), 2, fpEnvF);
	fbuf[0] = fineContainerEnd[0];
	fbuf[1] = fineContainerEnd[1];
	fwrite(fbuf, sizeof(float), 2, fpEnvF);
	for (int i = 0; i < ibuf[0]; i++)
	{
		fbuf[0] = fineWorld->GetBoundaryParticle(i)->m_curPosition[0];
		fbuf[1] = fineWorld->GetBoundaryParticle(i)->m_curPosition[1];
		fwrite(fbuf, sizeof(float), 2, fpEnvF);
	}

	ibuf[0] = coarseWorld->GetNumOfBoundaryParticles();
	fwrite(ibuf, sizeof(int), 1, fpEnvC);
	fbuf[0] = coarseContainerStart[0];
	fbuf[1] = coarseContainerStart[1];
	fwrite(fbuf, sizeof(float), 2, fpEnvC);
	fbuf[0] = coarseContainerEnd[0];
	fbuf[1] = coarseContainerEnd[1];
	fwrite(fbuf, sizeof(float), 2, fpEnvC);
	for (int i = 0; i < ibuf[0]; i++)
	{
		fbuf[0] = coarseWorld->GetBoundaryParticle(i)->m_curPosition[0];
		fbuf[1] = coarseWorld->GetBoundaryParticle(i)->m_curPosition[1];
		fwrite(fbuf, sizeof(float), 2, fpEnvC);
	}

	fclose(fpEnvC);
	fclose(fpEnvF);	
	printf("EnvForCoarse & FIne are saved.\n");
}

void CreateCoarseBreakingDam(std::vector<Vector2f>& p_damParticles)
{
	p_damParticles.resize(coarseDamWidth*coarseDamHeight);

	float diam = 2.0f * coarseR;
	float startX = -0.5f * containerWidth + diam + diam;
	float startY = diam + diam + diam;
	float yshift = sqrt(3.0f) * coarseR;

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int j = 0; j < coarseDamHeight; j++)
		{
			for (int i = 0; i < coarseDamWidth; i++)
			{
				p_damParticles[j*coarseDamWidth + i] = diam*Eigen::Vector2f((float)i, (float)j) + Eigen::Vector2f(startX, startY);
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
	p_damParticles.resize(fineDamWidth*fineDamHeight);

	float diam = 2.0f * fineR;
	float coarseDiam = 2.0f * coarseR;
	float startX = -0.5f * containerWidth + coarseDiam + coarseDiam;
	float startY = coarseDiam + coarseDiam + coarseDiam;
	float yshift = sqrt(3.0f) * fineR;

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int j = 0; j < fineDamHeight; j++)
		{
			for (int i = 0; i < fineDamWidth; i++)
			{
				p_damParticles[j*fineDamWidth + i] = diam*Eigen::Vector2f((float)i, (float)j) + Eigen::Vector2f(startX, startY);
			}
		}
	}
}
void CreateFineContainer(std::vector<Vector2f>& p_boundaryParticles)
{
	float x1 = -containerWidth / 2.0f;
	float x2 = containerWidth / 2.0f;
	float y1 = 0.0f;
	float y2 = containerHeight;

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

void ParticleTrainingDataSave()
{
	std::vector<FParticle2D*>& fineP = fineWorld->GetParticleList();
	std::vector<FParticle2D*>& coarseP = coarseWorld->GetParticleList();

	string cNeiInfoPath = "./PBFC2D_SD4/cNeiInfo_";
	string fNeiInfoPath = "./PBFC2D_SD4/fNeiInfo_";
	string GTPath = "./PBFC2D_SD4/GT_";
	string frameIdx = std::to_string(accFrameCount) + ".dat";
	
	FILE* fpCNei = fopen((cNeiInfoPath +frameIdx).c_str(), "wb");
	FILE* fpFNei = fopen((fNeiInfoPath + frameIdx).c_str(), "wb");
	FILE* fpGT = fopen((GTPath + frameIdx).c_str(), "wb");

	float fbuf[2];
	int ibuf[1];

	// num of fine particles 
	ibuf[0] = (int)fineP.size();
	fwrite(ibuf, sizeof(int), 1, fpCNei);
	fwrite(ibuf, sizeof(int), 1, fpGT);

	int maxCountC, maxCountCB, maxCountF, maxCountFB;

	for (int i = 0; i < fineP.size(); i++)
	{
		maxCountC = maxCountCB = maxCountF = maxCountFB = 0;

		// fine's GT
		Vector2f GT = fineP[i]->m_velocity;
		fbuf[0] = GT[0]; fbuf[1] = GT[1];
		fwrite(fbuf, sizeof(float), 2, fpGT);

		// fine's old velocity
		Vector2f deltaP = fineP[i]->m_tempVelocity;
		fbuf[0] = deltaP[0]; fbuf[1] = deltaP[1];
		fwrite(fbuf, sizeof(float), 2, fpCNei);

		// neighbor coarse particle info
		ibuf[0] = pbfc2D->m_tDataForCoarse[i].size();
		maxCountC < ibuf[0] ? maxCountC = ibuf[0] : maxCountC = maxCountC;
		fwrite(ibuf, sizeof(int), 1, fpCNei);
		for (int j = 0; j < ibuf[0]; j++)
		{
			fbuf[0] = pbfc2D->m_tDataForCoarse[i][j].mass;
			fwrite(fbuf, sizeof(float), 1, fpCNei);

			Vector2f& RVec = pbfc2D->m_tDataForCoarse[i][j].RVec;
			fbuf[0] = RVec[0]; fbuf[1] = RVec[1];
			fwrite(fbuf, sizeof(float), 2, fpCNei);

			Vector2f& dPos = pbfc2D->m_tDataForCoarse[i][j].dPos;
			fbuf[0] = dPos[0]; fbuf[1] = dPos[1];
			fwrite(fbuf, sizeof(float), 2, fpCNei);
		}

		// neighbor fine particle info
		ibuf[0] = pbfc2D->m_tDataForFine[i].size();
		maxCountF < ibuf[0] ? maxCountF = ibuf[0] : maxCountF = maxCountF;
		fwrite(ibuf, sizeof(int), 1, fpFNei);
		for (int j = 0; j < ibuf[0]; j++)
		{
			fbuf[0] = pbfc2D->m_tDataForFine[i][j].mass;
			fwrite(fbuf, sizeof(float), 1, fpFNei);

			Vector2f& RVec = pbfc2D->m_tDataForFine[i][j].RVec;
			fbuf[0] = RVec[0]; fbuf[1] = RVec[1]; 
			fwrite(fbuf, sizeof(float), 2, fpFNei);

			Vector2f& dPos = pbfc2D->m_tDataForFine[i][j].dPos;
			fbuf[0] = dPos[0]; fbuf[1] = dPos[1];
			fwrite(fbuf, sizeof(float), 2, fpFNei);
		}
	}
	
	//maxCountC > maxmaxCountCF ? maxmaxCountCF = maxCountC : maxmaxCountCF = maxmaxCountCF;
	//maxCountCB > maxmaxCountCB ? maxmaxCountCB = maxCountCB : maxmaxCountCB = maxmaxCountCB;
	//maxCountF > maxmaxCountFF ? maxmaxCountFF = maxCountF : maxmaxCountFF = maxmaxCountFF;
	//maxCountFB > maxmaxCountFB ? maxmaxCountFB = maxCountFB : maxmaxCountFB = maxmaxCountFB;

	//printf("maxmaxCF: %d, maxmaxCB: %d, maxmaxFF: %d, maxmaxFB: %d\n", maxmaxCountCF, maxmaxCountCB, maxmaxCountFF, maxmaxCountFB);
	
	fclose(fpCNei);
	fclose(fpFNei);
	fclose(fpGT);
}

void BoundaryParticleInfoSave()
{
	float fbuf[2];
	int ibuf[1];

	FILE* fpEnvF = fopen("./PBFC2D_SD10/BoundaryEnv.dat", "wb");

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
void ParticleTrainingDataSave_fineOnly_simulBefore()
{
	std::vector<FParticle2D*>& fineP = fineWorld->GetParticleList();
	FluidKernel2D& fineKernel = fineWorld->GetKernel();

	string frameIdx = std::to_string(accFrameCount) + ".dat";

	string NeighborIdxPath = "./PBFC2D_SD10/NeighborIdx_";
	string FParticleInfoPath = "./PBFC2D_SD10/FParticleInfo_";

	FILE* fpNeiIdx  = fopen((NeighborIdxPath + frameIdx).c_str(), "wb");
	FILE* fpFParticle = fopen((FParticleInfoPath + frameIdx).c_str(), "wb");
	
	float fbuf[2];
	int ibuf[1];

	// num of fine particles 
	ibuf[0] = (int)fineP.size();
	fwrite(ibuf, sizeof(int), 1, fpFParticle);

	int nNei, maxCountNei = 0;
		
	for (int i = 0; i < (int)fineP.size(); i++)
	{	
		FParticle2D* pi = fineP[i];
		
		// particle information save (mass, position, velocity)
		fbuf[0] = pi->m_mass;
		fwrite(fbuf, sizeof(float), 1, fpFParticle);
		fbuf[0] = pi->m_tempPosition[0];	fbuf[1] = pi->m_tempPosition[1];
		fwrite(fbuf, sizeof(float), 2, fpFParticle);
		fbuf[0] = pi->m_tempVelocity[0];	fbuf[1] = pi->m_tempVelocity[1];
		fwrite(fbuf, sizeof(float), 2, fpFParticle);
		
		// neighbor particle info (neighborList Size, neighbor Pid, neighbgor Pidx)
		ibuf[0] = nNei = (int)pi->m_neighborList.size();
		fwrite(ibuf, sizeof(int), 1, fpNeiIdx);
		for (int j = 0; j < (int)pi->m_neighborList.size(); j++)
		{
			if (pi->m_neighborList[j]->m_pid == Pid::Fluid)
			{
				ibuf[0] = 1;
				fwrite(ibuf, sizeof(int), 1, fpNeiIdx);
				ibuf[0] = pi->m_neighborList[j]->m_pIdx;
				fwrite(ibuf, sizeof(int), 1, fpNeiIdx);
			}
			else if (pi->m_neighborList[j]->m_pid == Pid::Boundary)
			{
				ibuf[0] = 0;
				fwrite(ibuf, sizeof(int), 1, fpNeiIdx);
				ibuf[0] = pi->m_neighborList[j]->m_pIdx;
				fwrite(ibuf, sizeof(int), 1, fpNeiIdx);
			}
		}
				
		maxCountNei < nNei ? maxCountNei = nNei : maxCountNei = maxCountNei;
	}
	
	maxmaxFB < maxCountNei ? maxmaxFB = maxCountNei : maxmaxFB = maxmaxFB;
	
	printf("(%d)frame maxmaxFB: %d \n", accFrameCount, maxmaxFB);

	fclose(fpNeiIdx);
	fclose(fpFParticle);
}
void ParticleTrainingDataSave_fineOnly_simulAfter()
{
	std::vector<FParticle2D*>& fineP = fineWorld->GetParticleList();
	FluidKernel2D& fineKernel = fineWorld->GetKernel();

	string GTVelPath = "./PBFC2D_SD10/GTVel_";
	string GTPosPath = "./PBFC2D_SD10/GTPos_";
	string frameIdx = std::to_string(accFrameCount) + ".dat";

	FILE* fpGTVel = fopen((GTVelPath + frameIdx).c_str(), "wb");
	FILE* fpGTPos = fopen((GTPosPath + frameIdx).c_str(), "wb");
	
	float fbuf[2];
	int ibuf[1];

	// num of fine particles 
	ibuf[0] = (int)fineP.size();
	fwrite(ibuf, sizeof(int), 1, fpGTVel);
	fwrite(ibuf, sizeof(int), 1, fpGTPos);
	
	for (int i = 0; i < (int)fineP.size(); i++)
	{
		FParticle2D* pi = fineP[i];
		Vector2f dPosFromF = fineWorld->GetPBFWorld2D()->m_deltaFromF[i];
		fbuf[0] = dPosFromF[0]; fbuf[1] = dPosFromF[1];
		fwrite(fbuf, sizeof(float), 2, fpGTPos);
		
		Vector2f corrVelFromF = dPosFromF / fineWorld->GetTimeStep();
		fbuf[0] = corrVelFromF[0]; fbuf[1] = corrVelFromF[1];
		fwrite(fbuf, sizeof(float), 2, fpGTVel);

		Vector2f dPosFromB = fineWorld->GetPBFWorld2D()->m_deltaFromB[i];
		fbuf[0] = dPosFromB[0]; fbuf[1] = dPosFromB[1];
		fwrite(fbuf, sizeof(float), 2, fpGTPos);

		Vector2f corrVelFromB = dPosFromB / fineWorld->GetTimeStep();
		fbuf[0] = corrVelFromB[0]; fbuf[1] = corrVelFromB[1];
		fwrite(fbuf, sizeof(float), 2, fpGTVel);
	}
	fclose(fpGTVel);
	fclose(fpGTPos);
}
*/