#pragma once

#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "gnuplot-iostream.h"

#include "signal_statistics.h"
#include "convolution.h"
#include "dft.h"
#include "idft.h"
#include "polar_notation.h"
#include "complex_dft.h"
#include "low_pass_windowed_sinc_filter.h"
#include "high_pass_windowed_sinc_filter.h"
#include "band_pass_windowed_sinc_filter.h"

//GLOBAL DEFINITIONS
//#define XP_SIG_LENGTH 4716 // test1;
#define XP_SIG_LENGTH 5186
#define SIG_LENGTH 320
#define IMP_RSP_LENGTH 29
#define ECG_SIG_LENGTH 640
#define POL_SIG_LENGTH 320
#define COMPLEX_SIG_LENGTH 501
#define KERNEL_LENGTH 320 

constexpr double SP7_ZERO_POSE[6] = {0.0, 0.0, 0.401, 0.0, 0.0, 0.0};
//SIGNALS in WAVEFORMS.CPP
extern double xplane_ax[XP_SIG_LENGTH];
//test 3
extern double xplane_t_test3[XP_SIG_LENGTH];
extern double xplane_ax_test3[XP_SIG_LENGTH];
extern double xplane_ay_test3[XP_SIG_LENGTH];
extern double xplane_az_test3[XP_SIG_LENGTH];
extern double xplane_vroll_test3[XP_SIG_LENGTH];
extern double xplane_vpitch_test3[XP_SIG_LENGTH];
extern double xplane_vyaw_test3[XP_SIG_LENGTH];

//
extern double InputSignal_f32_1kHz_15kHz[SIG_LENGTH];
extern double  Impulse_response[IMP_RSP_LENGTH];
extern double _640_points_ecg_[ECG_SIG_LENGTH];
extern double _320_pts_ecg_IMX[POL_SIG_LENGTH];
extern double _320_pts_ecg_REX[POL_SIG_LENGTH];
extern double _501pts_20Hz_sig_imx[COMPLEX_SIG_LENGTH];
extern  double _501pts_20Hz_sig_rex[COMPLEX_SIG_LENGTH];

//GLOBAL VARIABLE
double Output_signal_xp_rex[XP_SIG_LENGTH / 2];
double Output_signal_xp_imx[XP_SIG_LENGTH / 2];
double Output_signal_xp_mag[XP_SIG_LENGTH / 2];
double Output_signal_xp_pol_mag[XP_SIG_LENGTH / 2];
double Output_signal_xp_pol_phase[XP_SIG_LENGTH / 2];
double Output_signal_xp_dest[XP_SIG_LENGTH + KERNEL_LENGTH];
double Output_signal_xp_dest2[XP_SIG_LENGTH];
double Output_signal_xp_kernel[KERNEL_LENGTH];


double Output_signal2_xp_rex[XP_SIG_LENGTH / 2];
double Output_signal2_xp_imx[XP_SIG_LENGTH / 2];
double Output_signal2_xp_pol_mag[XP_SIG_LENGTH / 2];
double Output_signal2_xp_pol_phase[XP_SIG_LENGTH / 2];


double Output_signal_xp_integral[XP_SIG_LENGTH ];
double Output_signal_xp_double_integral[XP_SIG_LENGTH];

double MEAN;
double VARIANCE;
double STD;

double Output_signal[SIG_LENGTH + IMP_RSP_LENGTH];
double Output_signal_running_sum[SIG_LENGTH];
double Output_signal_dft_rex[SIG_LENGTH/2];
double Output_signal_dft_imx[SIG_LENGTH/2];
double Output_signal_dft_mag[SIG_LENGTH/2];
double Output_signal_idft[SIG_LENGTH];

double Output_signal_ecg_dft_rex[ECG_SIG_LENGTH / 2];
double Output_signal_ecg_dft_imx[ECG_SIG_LENGTH / 2];
double Output_signal_ecg_idft[ECG_SIG_LENGTH];

double Output_signal_pol_mag[POL_SIG_LENGTH];
double Output_signal_pol_phase[POL_SIG_LENGTH];

double Output_signal_cdft_freq_rex[COMPLEX_SIG_LENGTH];
double Output_signal_cdft_freq_imx[COMPLEX_SIG_LENGTH];

double Output_signal_lp_dest[SIG_LENGTH + KERNEL_LENGTH];
double Output_signal_lp_kernel[KERNEL_LENGTH];

//band pass filter
double Output_signal_bp_dest[SIG_LENGTH + KERNEL_LENGTH];
double Output_signal_bp_kernel[KERNEL_LENGTH];
double state_lw_cutoff_buff[KERNEL_LENGTH];
double state_up_cutoff_buff[KERNEL_LENGTH];




//FUNCTION PROTOTYPES IN MAIN FOR TESTING PURPOSE
void cueing_acceleration(double*, double, double, std::string, int);
void cueing_velocity(double*, double , double , std::string , int );
double Intergration_Trapezoidal(double input_curr, double input_prev, double output_prev, double t_prev, double t_curr);


//temp tests
void calc_running_sum(double*, double* , int );
void running_sum_test();
void dft_test();
void plot(std::vector<double> d1, std::vector<double> d2, std::vector<double> d3, std::vector<double> d4);
void idft_test();
void ecg_dft_idft_test();
void convolution_test();
void signal_statistics_test();
void r2p_test();
void complex_dft_test();
void low_pass_windowed_sinc_filter_test();

void band_pass_test();


//%% Boundary parameters
//
//% Home position-- > all the cranks are at the horizontal position
//zeroPose = [0 0 0.401 0 0 0];
//
//% Boundary translation limits
//
//xMinLimit = -0.05;
//yMinLimit = -0.05;
//zMinLimit = 0.396;
//xMaxLimit = 0.05;
//yMaxLimit = 0.05;
//zMaxLimit = 0.451;
//
//% Boundary rotational limits
//thetaxMinLimit = degtorad(-5);
//thetayMinLimit = degtorad(-5);
//thetazMinLimit = degtorad(-180);
//thetaxMaxLimit = degtorad(5);
//thetayMaxLimit = degtorad(5);
//thetazMaxLimit = degtorad(180);