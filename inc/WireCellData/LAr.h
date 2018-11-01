
namespace WireCell{

  class LAr{
  public: 
    LAr();
    
    // density of liquid Argon vs. temperature g/cm^3
    double Ldensity(double T);
    // density of gas Argon vs. temperature
    double Gdensity(double T);
    // vapor pressure vs. temperature bar
    double VPressure(double T);
    //Enthalpy: measure of the total energy of a thermodynamics system
    // enthalpy of vaporization in kJ/mol
    double Enthalpy(double T); 
    // Heat Capacity Cp, Cv
    double Cp(double T);
    double Cv(double T);
    //Speed of Sound 
    double SpeedofSound(double T);
    //Isothermal Compressibility  cm^2/dyne  1 Newton = 10^5 dyne
    // -1/V * dV/dp |_T
    double IsoCom(double T);
    // Viscosity  muP *s --> Force = Viscosity * A (area) * velocity / y (separation distance)
    double Viscosity(double T);
    
    //electron lifetime in ms, give a concentration in ppb
    double ele_lifetime(double con_ppb, double T=89, double E=0.5, int flag=1);
    double ks_f1(double x, double p0, double p1, double p2, double q1);
    double ks_f2(double x, double p0, double p1, double p2, double q1, double q2);
    double ks_f3(double x, double p0, double p1, double p2, double q1, double q2, double q3);
    
    // Boiling Line
    double BoilAr(double T);
    // Melting Line
    double MeltAr(double T);
    // Sublimation Line
    double SubLimeAr(double T);
    
    
    // Velocity of Ion drift  (theoretical)
    double vDarIonV(double T, double E=0.5);
    // Velocity of Ion drift (best fit) m/s
    double vDarIon(double T, double E=0.5);
    //Dielectric constant
    double epsilon(double T);
    //Index of refraction in Gas
    double GInrf(double lambda);
    //Index of refraction in Liquid
    double LInrf(double lambda, double T);
    // Rayleigh Scattering length cm
    double RRLAr(double lambda, double T);
    // Drift velocity
    double vD(double T, double E=0.5, int flag=1);
    // Drift velocity set I  unit: mm/us
    double vDrift(double T, double E=0.5);
    //Diffusions flag==1 logintudinal ==2 transverse
    // the sigma_L = sqrt(2*results*Delta_z / E), unit is eV
    double Diffusion(double T, double E=0.5, int flag=1);
    double fPade(double T, double t0, double p1, double p2, double q1, double q2);
    double fPoly(double T, double t0, double p1, double p2, double p3);
    //Recombinations Birks
    double recombine_Birks(double dEodx, double T, double E=0.5, int flag=1);
    //Recombinations Box model
    double recombine_Box(double dEodx, double T, double E=0.5, int flag=1);
    //Electron attachments
    
    
    void print_critical();
    void print_triple();
    void print_boiling();
    
    
    double T_c, p_c, rho_c; //Critical point
    // Critical point occurs under conditions at which no phase boundaries exist
    
    double T_t, p_t; // Triple Point
    // Three phases coexist in equilibruim
    
    double T_NBP;  // Normal Boiling Point
  };

};
