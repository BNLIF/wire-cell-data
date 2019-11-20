#include "WCPData/LAr.h"

#include "TMath.h"
#include <iostream>

using namespace WCP;
using namespace std;


double WCP::LAr::ele_lifetime(double con_ppb, double T, double E, int flag){
  double density = Ldensity(T)/units::g*pow(units::cm,3);
  double Mol = 1e-9 * 1000./39.948*density*con_ppb; // Mole is mole/L
 
  double results;
  double p1, p2, p0, q1,q2,q3;
  if (flag==1){
    //O2
    p0 = 10.5435;
    p1 = 9701.70895;
    p2 = 6020.07535;
    q1 = 958.28688;
    q2 = 415.2558;
    q3 = 1053.92648;
    results = 1e10*ks_f3(E/units::kilovolt*units::cm,p0,p1,p2,q1,q2,q3);
  }else if (flag==2){
    //CO2
    p0 = 0.529678;
    p1 = 0.504343;
    p2 = 4.481952;
    q1 = 1.69546;
    results = 1e10*ks_f1(E/units::kilovolt*units::cm,p0,p1,p2,q1);
  }else if (flag==3){
    //N2O
    p0 = 1.09709;
    p1 = 6.68216;
    p2 = -0.05016;
    q1 = 0.0276336;
    results = 1e11 * ks_f1(E/units::kilovolt*units::cm,p0,p1,p2,q1);
  }else if (flag==4){
    //H2O
    p0 = 1.09709;
    p1 = 6.68216;
    p2 = -0.05016;
    q1 = 0.0276336;
    results = 1.97023*1e11 * ks_f1(E/units::kilovolt*units::cm,p0,p1,p2,q1);
  }else if (flag==5){
    //SF6
    p0 = 1.884705;
    p1 = 0.026195;
    p2 = 0.00009571;
    q1 = 0.240995;
    results = 1e14*ks_f1(E/units::kilovolt*units::cm,p0,p1,p2,q1);
  }else if (flag==6){
    //N2
    results = 2.66347e7*(1+E/units::kilovolt*units::cm)/(E/units::kilovolt*units::cm)/density*vDrift(T,E)/units::mm*units::microsecond;
  }

  results = 1./(results * Mol)*1000.;

  return results*units::millisecond;
}

// good ...

double WCP::LAr::BoilAr(double T){
  double results;
  double a1 = -5.9410;
  double a2 = 1.35539;
  double a3 = -0.464976;
  double a4 = -1.5399;
  double phi = (1-T/T_c);
  results = p_c * exp(T_c/T*(a1*phi+a2*pow(phi,1.5)
			     +a3*pow(phi,2)+a4*pow(phi,4.5)));
  return results;
}
double WCP::LAr::MeltAr(double T){
  double results;
  double a1 = -7476.2665;
  double a2 = 9959.0613;
  results = p_t * (1+a1*(pow(T/T_t,1.05)-1)+a2*(pow(T/T_t,1.275)-1));
  return results;
}



double WCP::LAr::SubLimeAr(double T){
  double results;
  double a1 = -11.391604;
  double a2 = -0.39513;
  double phi = (1-T/T_t);
  results = p_t * exp(T_t/T*(a1*phi+a2*pow(phi,2.7)));
  return results;
}


double WCP::LAr::ks_f1(double x, double p0, double p1, double p2, double q1){
  double results;
  results = (p0+p1*x+p2*x*x)/(1+q1*x);
  return results;
}
double WCP::LAr::ks_f2(double x, double p0, double p1, double p2, double q1, double q2){
  double results;
  results = (p0+p1*x+p2*x*x)/(1+q1*x+q2*x*x);
  return results;
}
double WCP::LAr::ks_f3(double x, double p0, double p1, double p2, double q1, double q2, double q3){
  double results;
  results = (p0+p1*x+p2*x*x)/(1+q1*x+q2*x*x+q3*x*x*x);
  return results;
}



double WCP::LAr::recombine_Box(double dEodx,double T, double E, int flag){
  double density = Ldensity(T)/(units::g/pow(units::cm,3));
  double alpha = 0.93;
  double beta=0.212;

  
  E/=units::kilovolt/units::cm;
  dEodx/=units::MeV/units::cm;
  
  double results;
  results = log(alpha+beta*dEodx/E/density)/(beta*dEodx/E/density);
  if(flag==2){
    results = 1-0.803*results;
  }
  return results;
}



double WCP::LAr::recombine_Birks(double dEodx, double T, double E, int flag){
  double A = 0.8;
  double k = 0.0486;
  double alpha = 0.803;
  double results;
  E/=units::kilovolt/units::cm;
  dEodx/=units::MeV/units::cm;
  double density = Ldensity(T)/(units::g/pow(units::cm,3));
  
  results = A/(1. + k/E*dEodx/density);
  if (flag ==2){
    results = 1-alpha*results;
  }
  if (E<0.1 ||E>1.0) cout << "Warning: E = " << E << " kV/cm is out of range!" << endl;
  if (dEodx < 1.5 *density || dEodx > 30*density) cout << "Warning: dE/dx = " << dEodx/density << " MeV/cm is out of range!" << endl; 

  return results;
}



double WCP::LAr::Diffusion(double T, double E, int flag){
  double tT = T_t;
  double tC = T_c;

  double kB = 8.6173324e-5*units::eV/units::kelvin;

  double kT0T = kB * tT;
  double kT0C = kB * tC;

  double alpha = 2.2-1.5077*Ldensity(T)/(units::g/pow(units::cm,3));
  
  double results;
  if (flag==1){ //transverse
    results = (1-alpha)*fPoly(E/units::volt*units::cm/tT,kT0T/units::eV,0.0026289636,0.00003204097,-2.28720038e-7) + alpha * fPade(E/units::volt*units::cm/tC,kT0C/units::eV,0.110259,0.0020160892,0.0878054884,0.00020144744);
    results *= 40./56.;
  }else{ //longitduinal
    results = (1-alpha)*fPoly(E/units::volt*units::cm/tT,kT0T/units::eV,0.0012176581,6.322883559e-6,-2.79545739e-8) + alpha * fPade(E/units::volt*units::cm/tC,kT0C/units::eV,0.05057951,0.002007786,0.5480257,0.0055775997);
    results *= 16.5/20.95;
  }
  
  return results*units::eV;
}







double WCP::LAr::fPade(double x, double t0, double p1, double p2, double q1, double q2){
  double results ;
  results = (t0 + p1*x + p2*x*x)/(1+q1*x+q2*x*x);
  return results;
}

double WCP::LAr::fPoly(double x, double t0, double p1, double p2, double p3){
  double results = t0 + p1*x+p2*x*x+p3*x*x*x;
  return results;
}

double WCP::LAr::vDrift(double T, double E){
  
  //E/= units::kilovolt/units::cm;
  double xFit = 0.0938163 - 0.0052563 *(T/units::kelvin-87.302) - 0.000146981 *pow(T/units::kelvin-87.302,2);
  double muFit = 5.183987 + 0.01447761 * (T/units::kelvin-87.302) - 0.0034972*pow(T/units::kelvin-87.302,2)
    -0.0005162374*pow(T/units::kelvin-87.302,3);
  double results;

  
  if (E<xFit*units::kilovolt/units::cm){
    results = E*muFit  * units::mm/units::microsecond;
  }else if (E<=0.619*units::kilovolt/units::cm){
    results = vD(T,E,2);
  }else if (E>=0.699*units::kilovolt/units::cm){
    results = vD(T,E,1);
  }else{
    results = ((E/(units::kilovolt/units::cm)-0.619)/0.08*vD(T,E,1) + (0.699-E/(units::kilovolt/units::cm))/0.08*vD(T,E,2));
  }
  return results;
}



double WCP::LAr::vD(double T, double E, int flag){
  double p1, p2, p3, p4, p5, p6,t0;
  if (flag==1){
    //WalkowiakParameterSet
    p1 = -0.01481;
    p2 = -0.0075;
    p3 = 0.141;
    p4 = 12.4;
    p5 = 1.627;
    p6 = 0.317;
    t0 = 90.371;
  }else if (flag==2){
    //ICarusFit
    p1 = -0.04640231;
    p2 = 0.0171171;
    p3 = 1.881246;
    p4 = 0.9940772;
    p5 = 0.0117183;
    p6 = 4.202141;
    t0 = 105.7491;
  }else if (flag==3){
    // ICarus + Kalinin
    p1 = -0.0443247;
    p2 = 0.02063811;
    p3 = 1.975295;
    p4 = 1.1205820;
    p5 = 0.00594563;
    p6 = 5.807014;
    t0 = 101.05398861;
  }else if (flag==4){
    //ICarus + Kalinin + Walkowiak
    p1 = -0.0464023;
    p2 = 0.0171171;
    p3 = 1.8812463;
    p4 = 0.9940772;
    p5 = 0.0117183;
    p6 = 4.202141;
    t0 = 105.7491;
  }else if (flag==5){
    //ICarus + Kalini + Walkowiak + Aprile
    p1 = -0.05215565;
    p2 = 0.0228202;
    p3 = 1.9413135;
    p4 = -1.1397996;
    p5 = 0.0068283;
    p6 = 4.9838381;
    t0 = 99.28381;
  }
  double results;
  E /= units::kilovolt/units::cm;
  results = (1+p1*(T/units::kelvin-t0)) * (p3*E*log(1+fabs(p4)/E)+p5*pow(E,p6)) + p2 * (T/units::kelvin-t0);
  return results*units::mm/units::microsecond;
}


double WCP::LAr::RRLAr(double lambda, double T){
  
  double c = 29979245800.*units::cm/units::second;
  double kB = 1.38065e-16 *1e-7*units::joule/ units::kelvin;
  double alphaCM = 0.103276;
  double alphaLL = 0.104471;
  
  double omega = 2*3.1415926*c/lambda;
  double n = LInrf(lambda,T);
  double epsilon = 1 + 3*(-1.+n)*(1+n)*alphaCM/(alphaCM-n*n*alphaCM+(2+n*n)*alphaLL);

  double results = 1./(pow(omega,4)/6./3.1415926/pow(c,4)
			 *kB * T* IsoCom(T)*pow((-1+epsilon)*(2+epsilon)/3.,2));
  return results;
}

WCP::LAr::LAr(){
  T_c = 150.687* units::kelvin; // K 
  p_c = 48.630*units::bar; // bar
  rho_c = 0.5356*units::g/pow(units::cm,3); // g/cm^3
  
  T_t = 83.8058 * units::kelvin; // K
  p_t = 0.68891*units::bar; //bar = 0.1 MPa
  
  T_NBP = 87.3 * units::kelvin; // K
}

void WCP::LAr::print_critical(){
  cout << "Properties of the critical point: " << endl; 
  cout << "    Temperature: " << T_c/units::kelvin << " K" << endl;
  cout << "    Pressure:    " << p_c/units::bar << " bar or " << p_c/units::bar/10. << " MPa" << endl;
  cout << "    Density:     " << rho_c/ ( units::g/pow(units::cm,3)) << " g/cm^3" << endl; 
}

void WCP::LAr::print_triple(){
  cout << "Properties of triple point: " << endl;
  cout << "    Temperature: " << T_t/units::kelvin << " K" << endl;
  cout << "    Pressure:    " << p_t/units::bar << " bar or " << p_t/units::bar/10. << " MPa" << endl;  
}

void WCP::LAr::print_boiling(){
  cout << "Normal Boiling Point" << endl;
  cout << "    Temperature: " << T_NBP/units::kelvin << " K" << endl; 
}


double WCP::LAr::LInrf(double lambda, double T){
  //lambda/=units::nm;
  double nG = GInrf(lambda);
  double rhoG = 1.0034*0.0017840 ;
  double rhoL = Ldensity(T)/(units::g/pow(units::cm,3));
  //std::cout << nG << " " << rhoG << " " << rhoL << std::endl;
  double results = sqrt((2+nG*nG)*rhoG + 2*(-1+nG*nG)*rhoL)/
    sqrt((2+nG*nG)*rhoG+rhoL-nG*nG*rhoL);

  if (lambda < 110*units::nm) {
    cout << "Warning: wavelength " << lambda/units::nm << " nm is out of range!" << endl; 
    return 0;
  }
  return results;
}



double WCP::LAr::GInrf(double lambda){
  lambda/=units::nm;
  double c0 = 1.2055e-2;
  double a1 = 0.2075;
  double b1 = 91.012;
  double a2 = 0.0415;
  double b2 = 87.892;
  double a3 = 4.3330;
  double b3 = 214.02;
  lambda = lambda / 1000.;
  double results = 1 + c0 *(a1/(b1-1./lambda/lambda) 
			      + a2/(b2-1./lambda/lambda)
			      + a3/(b3-1./lambda/lambda));
  return results;
}





double WCP::LAr::epsilon(double T){
  double a = 4.12568;
  double rho = Ldensity(T)/39.948/units::g*pow(units::cm,3);
  double results = (-1.-2.*a*rho)/(-1+a*rho);
  return results;
}


double WCP::LAr::Enthalpy(double T){
  T/=units::kelvin;
  double a = 7.98304;
  double b = -0.0481275;
  double c = -0.0047259;
  
  double H = (a+b*T)/(1+c*T);
  if (T>100 || T < 83.8) {
    cout << "Warning: Temperature " << T << " K is out of range" << endl; 
    return 0;
  }
  return H * 1000*units::joule/units::mole;
}

double WCP::LAr::IsoCom(double T){
  double rcp = Cp(T); // kJ/kg/K
  double rcv = Cv(T); // kJ/kg/K
  double rss = SpeedofSound(T); // m/s
  double rhoL = Ldensity(T); // g/cm^3

  double results = rcp/rcv/(rhoL*rss*rss);
  return results;
  //units cm^2 * dyne ...
}

double WCP::LAr::Cv(double T){
  T/=units::kelvin;
  double a = 0.29934286;
  double b = 8.953781;
  double c = -21.8602116;
  double results = (a+b/T)/(1.0+c/T);
  return results* 1000 * units::joule/ units::kilogram /units::kelvin; 
}

double WCP::LAr::Cp(double T){
  T/=units::kelvin;
  double a = 2.126910;
  double b = -0.02416936;
  double c = 0.00014438771;
  double results = a + b *T + c*T*T;
  return results * 1000 * units::joule/ units::kilogram /units::kelvin;
}

double WCP::LAr::SpeedofSound(double T){
  T/=units::kelvin;
  double a = 1327.6370386;
  double b = -7.4603194;
  double c = -0.002213558263;
  double results = (a+b*T)/(1.0+c*T);
  return results * units::m/units::second;
}

double WCP::LAr::Viscosity(double T){
  double a = -0.390176214;
  double b = -65.2768756;
  double c = 1.215505389;
  double d = 128.5931277;
  double e = -0.618990054;
  double t = T/83.8058/units::kelvin;
  double results = (b/t+d/t/t)/(a+c/t+e/t/t);
  return results * (1e-6 * units::pascal * units::second);
}


double WCP::LAr::vDarIonV(double T, double E){
  
  E = E * 1000. / (units::volt/units::cm);
  double results = 0.432 * E * 0.01 / (Viscosity(T)/(1e-6 * units::pascal * units::second));
  return results * units::mm/units::second;
}
double WCP::LAr::vDarIon(double T, double E){
  T = T/units::kelvin;
  E = E * 1000. / (units::volt/units::cm); // V/cm ...
  double a = -1.6024763e-3;
  double b = 1.157022e-5;
  double c = 2.86982898e-7;
  double results = E * 0.01 * (a + b*T + c*T*T);
  return results * units::mm/units::second;
}




double WCP::LAr::VPressure(double T){
  double a = -5.9409785;
  double b = 1.3553888;
  double c = -0.46497607;
  double d = -1.5399043;
  double t = (1-T/T_c);
  
  double p = p_c * exp(T_c/T *(a*t
				 + b*pow(t,1.5)
				 + c*pow(t,2)
				 + d*pow(t,4.5)
				 ));
  return p;
}

double WCP::LAr::Ldensity(double T){
  double a = 1.5004262;
  double b = -0.31381290;
  double c = 0.086461622;
  double d = -0.041477525;

  double t = (1-T/T_c);
  
  double rho = log(rho_c) + a * pow(t,0.334)  
    + b * pow(t,2./3.) 
    + c * pow(t,7./3.)
    + d * pow(t,4.);
  rho = exp(rho);
  return rho ;
}

double WCP::LAr::Gdensity(double T){
  double a = -1.70695656;
  double b = -4.02739448;
  double c = 1.55177558;
  double d = -2.30683228;
  
  double t = (1-T/T_c);
  double rho = log(rho_c) + a * pow(t,0.345)
    + b * pow(t,5/6.)
    + c * pow(t,1)
    + d * pow(t,13./3.);
  rho = exp(rho);
  return rho ;
  
}


