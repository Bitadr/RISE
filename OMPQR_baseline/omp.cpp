#include "omp.h"

/*
 * X is the input signal
 * D is the Dictionary
 * l is the number of valid atoms in D, max of l is L_MAX
 * k is the number of atom needed to be selected
 * epsilon is the error threshold of OMP
 * V is the output coeff
 * supp is the index of coeff
 * supp_len is the number of non-zero coeff
 */
void omp(float Qt[L_MAX][M], float X[M], int l, int k, float epsilon, float V[L_MAX], int supp[L_MAX], int &supp_len)
{

	float R[L_MAX][L_MAX];
	float q_x[L_MAX];
	bool selected[L_MAX];
	float res[M];
	float new_col[M];
	float norm_x = 0;

	compute_normx_mul:
	for(int i=0;i<M;i++)
	{
		norm_x += X[i]*X[i];
	}
	load_res1:
	for(int i=0;i<M;i++)
	{
		res[i] = X[i];
	}

	load_selected:
	for(int i=0;i<l;i++)
	{
#pragma HLS LOOP_TRIPCOUNT min=1 max=128
		selected[i] = false;
	}


	int n;

	main_loop1:
	for(n=0;n<k;n++)
	{
#pragma HLS LOOP_TRIPCOUNT min=128 max=128
		float max = 0;
		int idx = 0;
		D_res_1:
		for(int k=0;k<l;k++)
		{
#pragma HLS LOOP_TRIPCOUNT min=128 max=128
			if(!selected[k])
			{
				float dot_D_res = 0;
				D_res_1_1:
				for(int j=0;j<M;j++)
				{
					dot_D_res += Qt[k][j]*res[j];
				}

				if(FABS(dot_D_res) > max)
				{
					idx = k;
					max = FABS(dot_D_res);
				}
			}
		}

		if(max < (float)1E-10)
		{
			break;
		}

		supp[n] = idx;
		selected[idx] = true;

		load_new_col:
		for(int i=0;i<M;i++)
		{
			new_col[i] = Qt[idx][i];
		}




		fill_R1:
		for(int i=0;i<n;i++)
		{
#pragma HLS LOOP_TRIPCOUNT min=1 max=128 avg=64

			float dot_Q_new_col = 0;
			fill_R1_1:
			for(int j=0;j<M;j++)
			{
				dot_Q_new_col += Qt[supp[i]][j]*new_col[j];
			}

			fill_R1_2:
			for(int j=0;j<M;j++)
			{
				new_col[j] -= dot_Q_new_col * Qt[supp[i]][j];
			}
			R[i][n] = dot_Q_new_col;
		}

		float norm_new_col = 0;
		norm_new_col_1:
		for(int j=0;j<M;j++)
		{
			norm_new_col += new_col[j]*new_col[j];
		}

		norm_new_col = sqrt(norm_new_col);

		R[n][n] = norm_new_col;

		fill_Q1:
		for(int j=0;j<M;j++)
		{
			Qt[idx][j] = new_col[j]/norm_new_col;
		}

		float dot_q_res = 0;

		dot_Q_res1:
		for(int j=0;j<M;j++)
		{
			 dot_q_res += Qt[idx][j] * res[j];
		}


		update_res1:
		for(int j=0;j<M;j++)
		{
			res[j] -= Qt[idx][j] * dot_q_res;
		}


		float norm_res = 0;

		norm_res1:
		for(int j=0;j<M;j++)
		{
			norm_res += res[j]*res[j];
		}


		float norm_diff = norm_res / norm_x;
		if( norm_diff < epsilon)
		{
			break;
		}

	}


	supp_len = n;
	dot_Q_X1:
	for(int j=0;j<n;j++)
	{
#pragma HLS LOOP_TRIPCOUNT min=1 max=128 avg=64
		float dot_Q_X = 0;
		dot_Q_X1_1:
		for(int i=0;i<M;i++)
		{
			dot_Q_X += Qt[supp[j]][i] * X[i];
		}
		q_x[j] = dot_Q_X;
	}

	solve_R_V1:
	for(int i=n-1;i>=0;i--)
	{
#pragma HLS LOOP_TRIPCOUNT min=1 max=128 avg=64
		float dot_V_R = 0;
		solve_R_V1_1:
		for(int j=1;j<n-i;j++)
		{
#pragma HLS LOOP_TRIPCOUNT min=1 max=128 avg=64
			dot_V_R += R[i][i+j] * V[i+j];
		}
		V[i] = (q_x[i] - dot_V_R) /R[i][i];
	}
}
