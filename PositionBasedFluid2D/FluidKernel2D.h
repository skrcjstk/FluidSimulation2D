#pragma once
#ifndef __FLUID_KERNEL2D_H__
#define __FLUID_KERNEL2D_H__

#include <math.h>
#include <Eigen/Dense>

using namespace Eigen;

class FluidKernel2D
{
public:
	const float PIF = 3.141592f;
	float h;
	float h2;
	float m_k;
	float m_l;
	float res_0;

	void SetSmoothingRadius(float p_radius)
	{
		h = p_radius;
		h2 = h * h;
		m_k = 40.0f / (7.0f*PIF*h2);
		m_l = 240.0f / (7.0f*PIF*h2);
		res_0 = Cubic_Kernel(Vector2f(0.0f, 0.0f));
	}
	float Cubic_Kernel(Vector2f r)
	{
		float rLength = r.norm();
		float q = rLength / h;
		float res;

		if (q <= 0.5)
		{
			float q2 = q*q;
			float q3 = q2*q;

			res = m_k * (6.0f*q3 - 6.0f*q2 + 1.0f);
		}
		else
		{
			res = m_k * (2.0f*pow(1.0f - q, 3));
		}

		return res;
	}
	float Cubic_Kernel0()
	{
		return res_0;
	}
	Vector2f Cubic_Kernel_Gradient(Vector2f r)
	{
		Vector2f res;
		float rl = r.norm();
		float q = rl / h;

		if (rl > 1.0e-6)
		{
			Vector2f gradq = r * ((float) 1.0 / (rl * h));
			if (q <= 0.5f)
			{
				res = m_l*q*((float) 3.0*q - (float) 2.0)*gradq;
			}
			else
			{
				const float factor = 1.0f - q;
				res = m_l*(-factor*factor)*gradq;
			}
		}

		return res;
	}

	float GetSmoothingRadius() { return h; }

	/*
	float selfFac;
	float invH;

	void SetSmoothingRadius(float p_radius)
	{
		invH = 1.0f / p_radius;
		selfFac = 15.0f / (7.0f * PIF) * invH * invH;
		res_0 = Cubic_Kernel(Vector2f(0.0f, 0.0f));
	}

	float Cubic_Kernel(Vector2f r)
	{
		float rLength = r.norm();
		float k = rLength * invH;
		
		float res;
		if (0.0f <= k && k < 1.0f)
		{
			res = (2.0f / 3.0f) - k*k + 0.5f*k*k*k;
		}
		else if (1.0f <= k && k < 2.0f)
		{
			res = powf((2.0f - k), 3) / 6.0f;
		}
		else if( k >= 2.0f)
		{
			res = 0;
		}

		return res * selfFac;
	}
	float Cubic_Kernel0()
	{
		return res_0;
	}
	Vector2f Cubic_Kernel_Gradient(Vector2f r)
	{
		float rLength = r.norm();
		float k = rLength * invH;

		Vector2f res;

		if (rLength > 1.0e-6)
		{
			Vector2f gradq = r * ((float) 1.0 / (rLength * h)); 
			if (0.0f <= k && k < 1.0f)
			{
				res = gradq * (-2.0f*k + 1.5f*k*k);
			}
			else if (1.0f <= k && k < 2.0f)
			{
				res = gradq * (-0.5f*(2.0f-k)*(2.0f-k));
			}
			else if (k >= 2.0f)
			{
				res[0] = res[1] = 0.0f;
			}
		}

		return res * selfFac;
	}
	*/
};

#endif