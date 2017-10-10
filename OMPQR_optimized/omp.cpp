#include "omp.h"

/*
 * mode = 0 reduce temp0(0:M-1) and use temp1(0:M-1) as temporary and return result
 * mode = 1 add temp1(0:M-1) to temp0(M:2*M-1) res
 * mode = 2 add temp0(0:M-1) to temp1(M:2*M-1) new_col
 * mode = 3 reduce temp1(0:M-1) and use temp0(0:M-1) as temporary and return result
 */
float reduce_temp(float temp0[2* M], float temp1[2* M], int mode)
{
#pragma HLS INLINE
	if(mode==0 || mode==3)
	{
		int n = M/2;
		bool use_temp0 = (mode==0);
		float reduce_ret = 0;

		compute_normx_sum_1:
		while(n>0)
		{
			if(use_temp0)
			{
				compute_normx_sum_2_0:
				for(int i=0;i<n;i++)
				{
#pragma HLS LOOP_TRIPCOUNT min=1 max=128
#pragma HLS UNROLL factor=16
#pragma HLS PIPELINE
					temp1[i] = temp0[2*i] + temp0[2*i+1];
				}
			}
			else
			{
				compute_normx_sum_2_1:
				for(int i=0;i<n;i++)
				{
#pragma HLS LOOP_TRIPCOUNT min=1 max=128
#pragma HLS UNROLL factor=16
#pragma HLS PIPELINE
					temp0[i] = temp1[2*i] + temp1[2*i+1];
				}
			}
			use_temp0 = !use_temp0;
			n = n>>1;
		}
		if(use_temp0)
		{
			reduce_ret = temp0[0];
		}
		else
		{
			reduce_ret = temp1[0];
		}
		return reduce_ret;
	}
	else if(mode==1)
	{
		equal_add_1:
		for(int i=0;i<M;i++)
		{
#pragma HLS DEPENDENCE variable=temp0 array inter false
#pragma HLS UNROLL factor=16
#pragma HLS PIPELINE
			temp0[i+M] += temp1[i];
		}
	}
	else if(mode==2)
	{

		equal_add_2:
		for(int i=0;i<M;i++)
		{
#pragma HLS DEPENDENCE variable=temp1 array inter false
#pragma HLS UNROLL factor=16
#pragma HLS PIPELINE
			temp1[i+M] += temp0[i];
		}
	}
	return 0;
}

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
#pragma HLS RESOURCE variable=X core=RAM_1P_BRAM
#pragma HLS ARRAY_PARTITION variable=Qt cyclic factor=16 dim=2
#pragma HLS ARRAY_MAP variable=X horizontal

#pragma HLS ARRAY_MAP variable=supp horizontal
#pragma HLS ARRAY_MAP variable=V horizontal

	float R[L_MAX][L_MAX];
	float q_x[L_MAX];
	bool selected[L_MAX];

 #pragma HLS ARRAY_MAP variable=R horizontal
 #pragma HLS ARRAY_MAP variable=selected horizontal
 #pragma HLS ARRAY_MAP variable=q_x horizontal

	float temp0[2*M];
	float temp1[2*M];

	#pragma HLS ARRAY_PARTITION variable=temp0 cyclic factor=16 dim=1
#pragma HLS ARRAY_PARTITION variable=temp1 cyclic factor=16 dim=1

	float norm_x = 0;

	compute_normx_mul:
	for(int i=0;i<M;i++)
	{
#pragma HLS UNROLL  factor=16
#pragma HLS PIPELINE
		temp0[i] = X[i]*X[i];
	}

	norm_x = reduce_temp(temp0, temp1, 0);
	supp_len = norm_x;



	load_res1:
	for(int i=0;i<M;i++)
	{
#pragma HLS UNROLL  factor=16
#pragma HLS PIPELINE
		temp0[i+M] = X[i];
	}


	load_selected:
	for(int i=0;i<l;i++)
	{

#pragma HLS LOOP_TRIPCOUNT min=128 max=128
#pragma HLS UNROLL  factor=16
#pragma HLS PIPELINE
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
#pragma HLS UNROLL  factor=16
#pragma HLS PIPELINE
					temp1[j] = Qt[k][j]*temp0[j+M];
				}

				dot_D_res = reduce_temp(temp0, temp1, 3);
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
#pragma HLS UNROLL  factor=16
			temp1[i+M] = Qt[idx][i];
		}




		fill_R1:
		for(int i=0;i<n;i++)
		{
#pragma HLS LOOP_TRIPCOUNT min=1 max=128 avg=64

			float dot_Q_new_col = 0;
			fill_R1_1:
			for(int j=0;j<M;j++)
			{
#pragma HLS UNROLL  factor=16
#pragma HLS PIPELINE
				temp0[j] = Qt[supp[i]][j]*temp1[j + M];
			}

			dot_Q_new_col = reduce_temp(temp0, temp1, 0);


			fill_R1_2:
			for(int j=0;j<M;j++)
			{
#pragma HLS UNROLL  factor=16
#pragma HLS PIPELINE
				temp0[j] = -1 * dot_Q_new_col * Qt[supp[i]][j];
			}
			reduce_temp(temp0, temp1, 2);


			R[i][n] = dot_Q_new_col;
		}

		float norm_new_col = 0;
		norm_new_col_1:
		for(int j=0;j<M;j++)
		{
#pragma HLS UNROLL  factor=16
#pragma HLS PIPELINE
			temp0[j] = temp1[j + M]*temp1[j + M];
		}

		norm_new_col = sqrt(reduce_temp(temp0, temp1, 0));

		R[n][n] = norm_new_col;

		fill_Q1:
		for(int j=0;j<M;j++)
		{
#pragma HLS UNROLL  factor=16
#pragma HLS PIPELINE
			Qt[idx][j] = temp1[j + M]/norm_new_col;
		}

		float dot_q_res = 0;

		dot_Q_res1:
		for(int j=0;j<M;j++)
		{
#pragma HLS UNROLL  factor=16
#pragma HLS PIPELINE
			 temp1[j] = Qt[idx][j] * temp0[j + M];
		}
		dot_q_res = reduce_temp(temp0, temp1, 3);


		update_res1:
		for(int j=0;j<M;j++)
		{
#pragma HLS UNROLL  factor=16
#pragma HLS PIPELINE
			temp0[j] = -1 * Qt[idx][j] * dot_q_res;
		}
		reduce_temp(temp0, temp1, 1);


		float norm_res = 0;

		norm_res1:
		for(int j=0;j<M;j++)
		{
#pragma HLS UNROLL  factor=16
#pragma HLS PIPELINE
			temp1[j] = temp0[j + M]*temp0[j + M];
		}
		norm_res = reduce_temp(temp0, temp1, 3);


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

		dot_Q_X1_1:
		for(int i=0;i<M;i++)
		{
#pragma HLS PIPELINE
			temp0[i] = Qt[supp[j]][i] * X[i];
		}
		q_x[j] = reduce_temp(temp0, temp1, 0);
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
#pragma HLS PIPELINE
			dot_V_R += R[i][i+j] * V[i+j];
		}
		V[i] = (q_x[i] - dot_V_R) /R[i][i];
	}


}
