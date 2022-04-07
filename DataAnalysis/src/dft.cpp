#include "dft.h"

void calc_sig_dft(double* sig_src_arr, double* sig_dest_rex_arr, double* sig_dest_imx_arr, int sig_length)
{
	int i, j, k;
	double PI = 3.14159265359;
	//Initialization
	for (j = 0; j < sig_length/2; j++)
	{
		sig_dest_rex_arr[j] = 0;
		sig_dest_imx_arr[j] = 0;
	}

	for (k = 0; k < sig_length / 2; k++)
	{
		for (i = 0; i < sig_length / 2; i++)
		{
			sig_dest_rex_arr[k] += sig_src_arr[i] * cos(((2 * PI * k * i) / sig_length));
			sig_dest_imx_arr[k] += (- sig_src_arr[i]) * sin(((2 * PI * k * i) / sig_length));
		}
	}
}

void get_dft_output_mag(double* sig_dest_mag_arr, double* sig_dest_rex_arr, double* sig_dest_imx_arr, int sig_length)
{
	int k;
	for  (k = 0; k < sig_length/2; k++)
	{
		sig_dest_mag_arr[k] = sqrt(pow(sig_dest_rex_arr[k],2) + pow(sig_dest_imx_arr[k], 2));
	}
}