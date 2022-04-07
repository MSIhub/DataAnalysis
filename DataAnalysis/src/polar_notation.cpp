#include "polar_notation.h"

void rect_to_polar_conversion(double*sig_src_rex_arr, double*sig_src_imx_arr,double*sig_out_mag_arr, double*sig_out_phase_arr, int sig_length )
{
	double PI = 3.14159265358979f;
	int k;

	for (k = 0; k < sig_length; k++)
	{
		sig_out_mag_arr[k] = sqrt(powf(sig_src_rex_arr[k],2) + powf(sig_src_imx_arr[k], 2));
		if (sig_src_rex_arr[k] == 0)
		{
			sig_src_rex_arr[k] = powf(10, -20);//researched value
			sig_out_phase_arr[k] = atan(sig_src_imx_arr[k] / sig_src_rex_arr[k]);
		}
		if (sig_src_rex_arr[k] < 0 && sig_src_imx_arr[k] < 0)
		{
			sig_out_phase_arr[k] -= PI; 
		}
		if (sig_src_rex_arr[k] < 0 && sig_src_imx_arr[k]>=0)
		{
			sig_out_phase_arr[k] += PI;
		}

	}

}