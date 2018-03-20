#include "PIC2D.h"

void PIC2D::Initialize(FluidWorld2D* p_world, Vector2f p_origin, Vector2f p_bSize, Vector2i p_nCount, float p_rho)
{
	int np = p_world->GetNumOfParticles();
	AssignResultF.resize(np);

	rho = p_rho;
	origin = p_origin;

	ni = p_nCount[0];
	nj = p_nCount[1];
	dx = p_bSize[0] / (float)ni;
	dy = p_bSize[1] / (float)nj;

	u = PICArray2d::Array2d<float>(ni + 1, nj, 0.0f);
	v = PICArray2d::Array2d<float>(ni, nj + 1, 0.0f);

	int cellCount = nj * ni;
	geo.resize(cellCount);
	cellsForF.resize(cellCount);
	cellsForB.resize(cellCount);
	cells_centor_pos.resize(cellCount);
	cells_centor_uv_coord.resize(cellCount);
	
	for (int j = 0; j < nj; j++)
		for (int i = 0; i < ni; i++)
		{
			int idx = j*ni + i;
			Vector2f uv_coord = Vector2f((i + 0.5f) * dx, (j + 0.5f) * dy);
			cells_centor_pos[idx] = origin + uv_coord;
			cells_centor_uv_coord[idx] = Vector2f(uv_coord[0] / dx, uv_coord[1] / dy);
			geo[idx] = A;
		}

	printf("PIC Grid(%d, %d)\n", ni, nj);
	printf("PIC dx(%.2f) dy(%.2f)\n", dx, dy);
}

void PIC2D::clean()
{
	u.~Array2d();
	v.~Array2d();

	geo.clear();
	cellsForF.clear();
	cellsForB.clear();
	cells_centor_pos.clear();
	cells_centor_uv_coord.clear();
	AssignResultF.clear();
}

void PIC2D::AssignBoundary(std::vector<FParticle2D*>& p_BpList)
{
	int nBp = (int)p_BpList.size();
	for (int n = 0; n < nBp; n++)
	{
		FParticle2D* p = p_BpList[n];

		int pi = (int)round((p->m_curPosition[0] - origin[0]) / dx);
		int pj = (int)round((p->m_curPosition[1] - origin[1]) / dy);
		int i = pi >= 0 && pi < ni ? pi : -1;
		int j = pj >= 0 && pj < nj ? pj : -1;

		int idx = j*ni + i;
		cellsForB[idx].push_back(p);
		geo[idx] = B;
	}
}

void PIC2D::AssignCells(FluidWorld2D* p_world)
{
	int np = p_world->GetNumOfParticles();

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)
		for (int j = 0; j < nj; ++j)
			for (int i = 0; i < ni; ++i)
			{
				int idx = j*ni + i;
				cellsForF[idx].clear();

				// marking as air
				if (geo[idx] == F)
					geo[idx] = A;
			}

#pragma omp for schedule(static)
		for (int n = 0; n < np; n++)
		{
			FParticle2D *p = p_world->GetParticle(n);

			int pi = (int)round((p->m_curPosition[0] - origin[0]) / dx);
			int pj = (int)round((p->m_curPosition[1] - origin[1]) / dy);
			int i = pi >= 0 && pi < ni ? pi : -1;
			int j = pj >= 0 && pj < nj ? pj : -1;

			AssignResultF[n][0] = i;
			AssignResultF[n][1] = j;
		}
	}

	for (int n = 0; n < np; n++)
	{
		FParticle2D *p = p_world->GetParticle(n);
		Vector2i& assign = AssignResultF[n];
		if (assign[0] != -1 && assign[1] != -1)
		{
			int idx = assign[1] * ni + assign[0];
			cellsForF[idx].push_back(p);

			// marking as fluid
			if (geo[idx] != B)
				geo[idx] = F;
		}
	}
}

void PIC2D::Map_P2G(FluidWorld2D* p_world)
{
	float radii = p_world->GetParticleRadius();

#pragma omp parallel default(shared)
	{
		// u-component of velocity
#pragma omp for schedule(static)
		for (int j = 0; j < nj; j++)
			for (int i = 0; i < ni + 1; i++)
			{
				std::vector<FParticle2D *> neighbors;
				GetNeigboringParticles_cell(i, j, -1, 0, -1, 1, neighbors);

				if ((int)neighbors.size() > 0)
				{
					Vector2f pos = Vector2f(i * dx, (j + 0.5f)*dy) + origin;

					float sum_weight = 0.0;
					float sum_u = 0.0;
					for (FParticle2D* p : neighbors)
					{
						float weight = 4.0f / 3.0f * M_PI * rho * radii * radii * radii * linear_kernel(p->m_curPosition - pos, dx);
						sum_u += weight * p->m_velocity[0]; // +p->m_c.col(0).dot(pos - p->m_curPosition);
						sum_weight += weight;
					}

					if (sum_weight != 0.0)
						u.set(i, j, sum_u / sum_weight);
					else
						u.set(i, j, 0.0);
				}
				else
					u.set(i, j, 0.0);
			}

		// v-component of velocity
#pragma omp for schedule(static)
		for (int j = 0; j < nj + 1; j++)
			for (int i = 0; i < ni; i++)
			{
				std::vector<FParticle2D *> neighbors;
				GetNeigboringParticles_cell(i, j, -1, 1, -1, 0, neighbors);

				if ((int)neighbors.size() > 0)
				{
					Vector2f pos = Vector2f((i + 0.5)*dx, j * dy) + origin;

					float sum_weight = 0.0;
					float sum_v = 0.0;
					for (FParticle2D* p : neighbors)
					{
						float weight = 4.0f / 3.0f * M_PI * rho * radii * radii * radii * linear_kernel(p->m_curPosition - pos, dy);
						sum_v += weight * p->m_velocity[1]; //+ p->m_c.col(1).dot(pos - p->m_curPosition);
						sum_weight += weight;
					}

					if (sum_weight != 0.0)
						v.set(i, j, sum_v / sum_weight);
					else
						v.set(i, j, 0.0);
				}
				else
					v.set(i, j, 0.0);
			}
	}
}

void PIC2D::GetNeigboringParticles_cell(int i, int j, int wl, int wh, int hl, int hh, std::vector<FParticle2D *>& res)
{
	for (int sj = j + hl; sj <= j + hh; sj++)
		for (int si = i + wl; si <= i + wh; si++)
		{
			if (si < 0 || si > ni - 1 || sj < 0 || sj > nj - 1)
				continue;
			res.insert(res.end(), cellsForB[(sj * ni) + si].begin(), cellsForB[(sj * ni) + si].end());
			res.insert(res.end(), cellsForF[(sj * ni) + si].begin(), cellsForF[(sj * ni) + si].end());
		}
}

Vector2f PIC2D::GetVelocity(int p_gridIdx)
{
	//Interpolate the velocity from the u and v grids
	Vector2f& uv_coord = cells_centor_uv_coord[p_gridIdx];
	Vector2f p0 = uv_coord - Vector2f(0, 0.5);
	Vector2f p1 = uv_coord - Vector2f(0.5, 0);
	float u_value = interpolate_value(p0, u);
	float v_value = interpolate_value(p1, v);

	return Vector2f(u_value, v_value);
}

void PIC2D::GetDescriptorAll(float result[], int desc_width)
{
	int d = desc_width;
	int halfCnt = (int)(desc_width / 2);
	int dataSize = 3 * (d*d);	// m,u,v * (desc_width * desc_width)
	int np = (int)AssignResultF.size();

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)
		// m,u,v
		for (int n = 0; n < np; n++)
		{
			Vector2i& ij = AssignResultF[n];
			int startidx = n * dataSize;

			for (int j = -halfCnt; j <= halfCnt; j++)
				for (int i = -halfCnt; i <= halfCnt; i++)
				{
					int idx = 3 * ((j + halfCnt)*(d)+(i + halfCnt));
					
					Vector2i neiGrid = ij + Vector2i(i, j);
					int neiIdx = neiGrid[1] * ni + neiGrid[0];
					Vector2f vel = GetVelocity(neiIdx);

					result[startidx + idx + 0] = geo[neiIdx];
					result[startidx + idx + 1] = vel[0];
					result[startidx + idx + 2] = vel[1];
				}
		}
	}
}
void PIC2D::GetDescriptor(Vector2i p_ij, float result[], int desc_width)
{
	int d = desc_width;
	int halfCnt = (int)(desc_width / 2);

	// m,u,v
	for (int j = -halfCnt; j <= halfCnt; j++)
		for (int i = -halfCnt; i <= halfCnt; i++)
		{
			int idx = 3 * ((j + halfCnt)*(d)+(i + halfCnt));

			Vector2i neiGrid = p_ij + Vector2i(i, j);
			int neiIdx = neiGrid[1] * ni + neiGrid[0];
			Vector2f vel = GetVelocity(neiIdx);

			result[idx + 0] = geo[neiIdx];
			result[idx + 1] = vel[0];
			result[idx + 2] = vel[1];
		}

}
inline float PIC2D::interpolate_value(Vector2f& point, PICArray2d::Array2d<float>& grid)
{
	int i, j;
	float fx, fy;
	float result;

	get_barycentric(point[0], i, fx, 0, grid.width);
	get_barycentric(point[1], j, fy, 0, grid.height);

	return bilerp(grid.get(i, j), grid.get(i + 1, j), grid.get(i, j + 1), grid.get(i + 1, j + 1), fx, fy);

}
inline Vector2f PIC2D::affine_interpolate_value(Vector2f& point, PICArray2d::Array2d<float>& grid)
{
	int i, j;
	float fx, fy;
	Vector2f result;

	get_barycentric(point[0], i, fx, 0, grid.width);
	get_barycentric(point[1], j, fy, 0, grid.height);

	result = grad_bilerp(grid.get(i, j), grid.get(i + 1, j), grid.get(i, j + 1), grid.get(i + 1, j + 1), fx, fy);

	return result;

}
void PIC2D::UpdateAffineMatrix(FluidWorld2D* p_world)
{
	for (int i = 0; i < p_world->GetNumOfParticles(); i++)
	{
		FParticle2D* particle = p_world->GetParticle(i);
		Vector2f dist = (particle->m_curPosition - origin);
		Vector2f p(dist[0] / dx, dist[1] / dy);
		Vector2f p0 = p - Vector2f(0, 0.5);
		Vector2f p1 = p - Vector2f(0.5, 0);

		particle->m_c.col(0) = affine_interpolate_value(p0, u) / dx;
		particle->m_c.col(1) = affine_interpolate_value(p1, v) / dy;
	}
}
