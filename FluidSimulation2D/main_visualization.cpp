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

void loadData();

Primitive spherePrimiCoarse, spherePrimiFine;
GLint context_major_version, context_minor_version;
bool doPause = true;

float fineR = 0.025f;
float coarseR = 0.05f;

int frameLimit = 3000;
int accFrameCount = 0;
float accTime = 0.0f;

std::vector<Vector2f> fineP;
std::vector<Vector2f> coarseP;

int main(int argc, char** argv)
{
	// OpenGL
	MiniGL::init(argc, argv, 1024, 768, 0, 0, "Fluid load");
	MiniGL::initLights();
	MiniGL::setClientIdleFunc(50, timeStep);
	MiniGL::setKeyFunc(0, 'r', reset);

	MiniGL::getOpenGLVersion(context_major_version, context_minor_version);

	MiniGL::setClientSceneFunc(render);
	MiniGL::setViewport(40.0f, 0.1f, 500.0f, Vector3f(0.0, 4.5, 13.0), Vector3f(0.0, 4.5, 0.0));

	TwAddVarRW(MiniGL::getTweakBar(), "Pause", TW_TYPE_BOOLCPP, &doPause, " label='Pause' group=Simulation key=SPACE ");

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

	// load data
	loadData();
	
	doPause != doPause;

	accFrameCount += 1;
	if (accFrameCount == frameLimit)
		doPause = true;
}

void reset()
{
}
void render()
{
	MiniGL::coordinateSystem();
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
		// drawing fine world
	{
		Vector3f translationForFine(-1.5f, 0.0f, 0.0f);
		float fluidColor[4] = { 0.8f, 0.3f, 0.3f, 1.0f };
		for (int i = 0; i < (int)fineP.size(); i++)
		{
			Vector3f particlePos = Vector3f(fineP[i][0], fineP[i][1], 0.0f) + translationForFine;
			spherePrimiFine.renderSphere(particlePos, fluidColor);
		}
	}
	}

	// drawing coarse world
	{
		Vector3f translationForCoarse(1.5f, 0.0f, 0.0f);
		float fluidColorForCoarse[4] = { 0.0f, 0.2f, 0.7f, 1.0f };
		for (int i = 0; i < (int)coarseP.size(); i++)
		{
			Vector3f particlePos = Vector3f(coarseP[i][0], coarseP[i][1], 0.0f) + translationForCoarse;
			spherePrimiCoarse.renderSphere(particlePos, fluidColorForCoarse);
		}
	}

	// drawing union world
	{
		Vector3f translationForUnion(0.0f, 1.5f, 0.0f);

		// drawing fine particles
		float fluidColor[4] = { 0.8f, 0.3f, 0.3f, 0.2f };
		for (int i = 0; i < (int)fineP.size(); i++)
		{
			Vector3f particlePos = Vector3f(fineP[i][0], fineP[i][1], 0.0f) + translationForUnion;
			spherePrimiFine.renderSphere(particlePos, fluidColor);
		}

		// drawing fine particles
		float fluidColorForCoarse[4] = { 0.0f, 0.2f, 0.7f, 0.5f };
		for (int i = 0; i < (int)coarseP.size(); i++)
		{
			Vector3f particlePos = Vector3f(coarseP[i][0], coarseP[i][1], 0.0f) + translationForUnion;
			spherePrimiCoarse.renderSphere(particlePos, fluidColorForCoarse);
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

void loadData()
{
	string positionDataPath = "E:/Tensorflow/tensorflow/tensorflow/contrib/cmake/build/modelLoader/PBFC2DTrainedModel/PBFC2DTrainResult_180528_training2/";
	string frameIdx = std::to_string(accFrameCount) + ".dat";

	string FinePositionPath = positionDataPath + "FinePos_";
	FILE* fpFineP = fopen((FinePositionPath + frameIdx).c_str(), "rb");

	string CoarsePositionPath = positionDataPath + "CoarsePos_";
	FILE* fpCoarseP = fopen((CoarsePositionPath + frameIdx).c_str(), "rb");

	int ibuf[1];
	float fbuf[2];

	fread(ibuf, sizeof(int), 1, fpFineP);
	int nFP = ibuf[0];
	fineP.resize(nFP);
	printf("%d nfp\n", nFP);
	for (int i = 0; i < nFP; i++)
	{
		fread(fbuf, sizeof(float), 2, fpFineP);
		fineP[i][0] = fbuf[0];
		fineP[i][1] = fbuf[1];
	}

	fread(ibuf, sizeof(int), 1, fpCoarseP);
	int nCP = ibuf[0];
	coarseP.resize(nCP);
	for (int i = 0; i < nCP; i++)
	{
		fread(fbuf, sizeof(float), 2, fpCoarseP);
		coarseP[i][0] = fbuf[0];
		coarseP[i][1] = fbuf[1];
	}
	printf("%d frame loaded.\n", accFrameCount);

	fclose(fpFineP);
	fclose(fpCoarseP);
}
*/