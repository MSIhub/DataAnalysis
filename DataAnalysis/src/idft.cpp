#include "idft.h"

void calc_sig_idft(double* sig_idft_out_arr, double* sig_src_rex_arr, double* sig_src_imx_arr, int idft_length)
{
    double PI = 3.14159265359;
    int i, k;
    for (k = 0; k < idft_length / 2; k++)
    {
        sig_src_rex_arr[k] = sig_src_rex_arr[k] / (idft_length/2);
        sig_src_imx_arr[k] = -sig_src_imx_arr[k] / (idft_length/2);
    }

    sig_src_rex_arr[0] = sig_src_rex_arr[0] / 2;
    sig_src_imx_arr[0] = -sig_src_imx_arr[0] / 2;

    for (i = 0; i < idft_length; i++)
    {
        sig_idft_out_arr[i] = 0;
    }

    for (k = 0; k < idft_length/2 ; k++)
    {
        for (i = 0; i < idft_length; i++)
        {
            auto real_part = sig_src_rex_arr[k] * cos(2 * PI * k * i / (double)idft_length);
            auto img_part =  sig_src_imx_arr[k] * sin(2 * PI * k * i / (double)idft_length);
            sig_idft_out_arr[i] += real_part + img_part;
            //sig_idft_out_arr[i] = sig_idft_out_arr[i] + sig_src_rex_arr[k] * cos(2 * PI * k * i / (double) idft_length);
            //sig_idft_out_arr[i] = sig_idft_out_arr[i] + sig_src_imx_arr[k] * sin(2 * PI * k * i / (double) idft_length);
        }
    }	
}