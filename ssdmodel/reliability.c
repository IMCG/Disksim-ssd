#include "reliability.h"
#include<stdio.h>
//physical constants.

#define pi (22/7)
static double electron_charge = 1.60217646e-19; //charge of electron
static double reduced_barrier_height = 1.1*1.60217646e-19; //charge of electron
static double barrier_height = 3.2*1.60217646e-19; //#3.2eV - in Joules 
static double plancks_constant = 6.626068e-34;// #planck's constant. - in J.s
static double electron_mass = 9.10938188e-31;// #mass of electron
static double reduced_electron_mass = 0.42 * 9.10938188e-31;// #reduced electron mass
static double permittivity_vacuum = 8.854187817e-12;
static double silicon_permitivity = 4.5 * 8.854187817e-12; //#permitivity of silicon
static double reduced_plancks_constant = 6.626068e-34/(2*3.1415926535897931);
static double boltzmann_constant = 1.380650424e-23;//# Boltzmann's constant. JK^-1


//device parameters.
static const double f = 80e-9;// #80nm - feature size.
static const double t_ox = 7.5e-9;//#thickness of oxide -7.5e-9
static const double t_ono = 12e-9;//#thickness of ono - 12nm.
static const double gcr = 0.6;//#gcr = c_ono/c_fg

//Device parameters used globally (in this file).
//initialized in initialize_device_paramters.
static double qcox = -1;
static double c_ono_area = -1;

//temperature parameters.
static double t_use = 278+25;
static double t_stress = 278+25;
static double silc_af = -1;
static double detrap_af = -1;
static double silc_ae = 0.3*1.60217646e-19;
static double detrap_ae  = 1.1*1.60217646e-19;

//voltage parameters
static double flat_band = 1;// #flat band voltage.
static double delvth_slc = 1.70;//#slc separation between levels.
static double delvth_mlc = 0.65;//#mlc separation between levels.
static double delvth = -1;


inline double arrhenius_factor(double t_use, double t_stress, double ae) 
{
	return exp((ae/boltzmann_constant)*(1.0/t_use-1.0/t_stress));
}

void initialize_device_physical_parameters()
{
	double wl = 1;
	double area = f * f * wl;
	double c_ox_area = silicon_permitivity/t_ox;
	double cox = c_ox_area * area;//#(e/d)*A - where a = f*f
	c_ono_area = (gcr*c_ox_area);
	qcox = electron_charge/cox;
	delvth = delvth_mlc;//defaults to MLC flash.

	silc_af = arrhenius_factor(t_use,t_stress,silc_ae);
	detrap_af = arrhenius_factor(t_use,t_stress,detrap_ae);
  fprintf(stderr,"DISKSIM SILC_af:%f\n",silc_af);
  fprintf(stderr,"DISKSIM detrap_af:%f\n",detrap_af);
  fprintf(stderr,"DISKSIM Delvth:%f\n",delvth);
  fprintf(stderr,"DISKSIM Qcox:%E\n",qcox);
  fprintf(stderr,"DISKSIM c_ono_area:%E\n",c_ono_area);

}

inline double trap(double cycles) {
	return (0.08 * pow(cycles,0.62) + 5 *pow(cycles,0.3))*qcox;
}

inline double trap_derivative(double cycles) {
	return (0.08*0.62*pow(cycles,0.62-1) + 5*0.3*pow(cycles,0.3-1))*qcox;
}

inline double detrap(double vth, double recovery_time) {
	return log(vth*1000)*log(recovery_time)/1000;
}


double find_root(double vth, double cycle)//the eqn cycle cannot be more than cycle.
{
  //Newton's method to find the root of a polynomial.
  double newcycle = cycle,prev_cycle=0;
  double threshold = 0.00001;
  double a=0;
  double dvth = trap(newcycle);
  while(fabs(newcycle-prev_cycle)>threshold && fabs(dvth-vth)>threshold)
  {
    prev_cycle = newcycle;
    a = (dvth-vth)/trap_derivative(newcycle);
    newcycle = (a>newcycle)?(a-newcycle):(newcycle-a);
    dvth = trap(newcycle);
  }
  if(fabs(trap(newcycle)-vth)<=threshold)
	  return newcycle;
	else
		return 0;
}

double calc_tunnel_current(double potential_diff, double trap_vth, double barrier_height)
{
  double f_ox = (potential_diff-trap_vth)/t_ox;
  double magic_no = ( pi * boltzmann_constant * t_use * 2 *
             sqrt(2*reduced_electron_mass*barrier_height) /
             (reduced_plancks_constant * electron_charge * f_ox) );// #check page 62 in the thin book.
  double sin_magic_no = sin(magic_no);
  double b_fn = ( (4* sqrt(2*reduced_electron_mass) *
            pow(barrier_height,1.5)) /
            (3*reduced_plancks_constant*electron_charge) );
  double a_fn = ( ((0.19/0.42) * pow(electron_charge,3)) /
          (8*pi * plancks_constant * barrier_height) );
  double c_fn = exp(-1*b_fn/f_ox);
  double tc = (magic_no/sin_magic_no) * a_fn * pow(f_ox,2) * c_fn;
  return tc;
}

double calc_retention_period(double trap_vth, double cycles)
{
	//constants - C,D,alpha and beta - Will give them better names later.
  double alpha = -0.4;
  double beta = 0.3;// #0.35
  double C =  4.61e-22;//#10.95e-22
  double D = 0.69;//#0.47
  double mean_vth = 2.9;// #2*delvth
  double pgm_voltage=16;
  double del_q_inj = mean_vth * c_ono_area;// #c_ox_area or cox?
  double q_inj = pow((del_q_inj * cycles),alpha);
  double j_stress = calc_tunnel_current(pgm_voltage*gcr,flat_band,barrier_height);

	double a_silc,b_silc,f_ox,j_ss,retention_period;
  //print "Stress current:%E(mA/cm2)" % (j_stress/1e1)
  //print "Q_inj:%E (C)" % (q_inj*area)

  if(cycles >=10 && cycles <= 100)
    reduced_barrier_height = 1.1*electron_charge;
  else if(cycles >100 && cycles<=1000) //#1.1eV-1.005eV
    reduced_barrier_height = (1.1-((cycles-101)*0.05/900))*electron_charge;
  else if(cycles >1000 && cycles<=10000) //#1.05-eV-1eV
    reduced_barrier_height = (1.05-((cycles-1001)*0.05/9000))*electron_charge;
  else if(cycles >10000 && cycles<=100000) //#1.00-eV-0.95eV
    reduced_barrier_height = (1.0-((cycles-10001)*0.05/90000))*electron_charge;
  else
    reduced_barrier_height = 0.9*electron_charge;

  //#reduced_barrier_height = 0.9*electron_charge #override all the above.    
  //#print reduced_barrier_height/electron_charge
  b_silc = ((4* sqrt(2*reduced_electron_mass) *
            pow(reduced_barrier_height,1.5)) /
            (3*reduced_plancks_constant*electron_charge));
  f_ox = mean_vth/t_ox;
  //print "Oxide field during retention:%E(MV/cm)" % (f_ox/1e8)
  //since unit of C is mA^(1-beta), convert j_stress to mA multiply by 1000)
  //apply beta to it and then multiply C and this result. Since the final result
  //is still in mA, divide by 1e3 to get the result in Amps.
  a_silc = ( C * pow(j_stress*1e3,beta) * exp(-1*D*q_inj) )/1e3;
  j_ss =  a_silc * pow(f_ox,2) * exp(-1*b_silc/f_ox);
  //print "j^beta:%E,qterm1:%E,qterm2:%E,A_silc:%E" % (x1,x2,x3,a_silc)
  //j_ss =  constant * a1 * a2 * d_ot * tunnel_current(dvth,reduced_barrier_height)
  //print "%E,%E" % (c_ox_area,c_ono_area)
  retention_period = ((delvth - trap_vth) * c_ono_area ) / j_ss;// #both c_ox_area and j_ss are per unit area
  retention_period = retention_period * silc_af; 
  //print "Leakage current:%E(A/cm^2),retention_period:%E" % (j_ss*1e-4,retention_period/a_year)
  return retention_period;
}
