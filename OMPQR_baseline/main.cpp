#include <iostream>
#include "omp.h"

using namespace std;

int main()
{
	int l = 128;

	float X[M];
	float Dt[L_MAX][M];
	float Qt[L_MAX][M];
	for(int j=0;j<l;j++)
	{
		float norm_d = 0;
		for(int i=0;i<M;i++)
		{
			Dt[j][i] = (rand()%256)/256.0;
			norm_d += Dt[j][i]*Dt[j][i];
		}
		for(int i=0;i<M;i++)
		{
			Dt[j][i] = Dt[j][i] / norm_d;
			Qt[j][i] = Dt[j][i];
		}
	}

	float norm_x = 0;
	for(int j=0;j<M;j++)
	{
		X[j] = (rand()%256)/256.0;
		norm_x += X[j]*X[j];
	}

	float V[L_MAX];
	int supp[L_MAX];
	int supp_len;


	omp(Qt, X, l, l, 0.13, V, supp, supp_len);

	float Xhat[M];
	for(int j=0;j<M;j++)
	{
		Xhat[j] = 0;
	}

	for(int i=0;i<supp_len;i++)
	{
		for(int j=0;j<M;j++)
		{
			Xhat[j] += Dt[supp[i]][j] * V[i];
		}
	}

	float norm_diff = 0;
	for(int j=0;j<M;j++)
	{
		norm_diff += (X[j] - Xhat[j])*(X[j] - Xhat[j]);
	}
	cout << " norm = " << norm_diff/ norm_x << endl;
	cout << " supp_len = " << supp_len << endl;
	cout << "v(i) i" << endl;
	for(int i=0;i<supp_len;i++)
	{
		cout << V[i] << " " << supp[i] << endl;
	}




	return 0;
}
