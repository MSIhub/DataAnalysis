#include "complex_dft.h"

void complex_dft(double *sig_src_time_domain_rex, double *sig_src_time_domain_imx, double*sig_dest_freq_domain_rex, double*sig_dest_freq_domain_imx, int sig_length)
{
	double PI = 3.14159265;
	double SR, SI;
	for (int k = 0; k < sig_length-1; k++)
	{
		for (int i = 0; i < sig_length-1; i++)
		{
			SR = cos((2 * PI * k * i) / sig_length);
			SI = -sin((2 * PI * k * i) / sig_length);
			
			sig_dest_freq_domain_rex[k] += (sig_src_time_domain_rex[i] * SR) - (sig_src_time_domain_imx[i] * SI);
			sig_dest_freq_domain_imx[k] += (sig_src_time_domain_imx[i] * SI) - (sig_src_time_domain_imx[i]*SR);
		}
	}
}