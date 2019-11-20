#ifndef WCPData_MCParticle_h
#define WCPData_MCParticle_h

namespace WCP{
  class MCParticle{
  public: 
    MCParticle(){};
    ~MCParticle(){};

    int pdg;
    float startXYZT[4];
    float endXYZT[4];
    float startMomentum[4];
    std::vector<Point> trajectory;
  };
  
  typedef std::vector<WCP::MCParticle*> MCParticleSelection;

}



#endif
