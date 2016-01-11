#ifndef WireCellData_MCParticle_h
#define WireCellData_MCParticle_h

namespace WireCell{
  class MCParticle{
    MCParticle(){};
    ~MCParticle(){};
  public: 
    int pdg;
    float startXYZT[4];
    float endXYZT[4];
    float startMomentum[4];
    std::vector<Point> trajectory;
  };
}


#endif
