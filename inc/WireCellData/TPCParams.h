namespace WireCell{
  class TPCParams {
    double m_pitch_u; // wire pitch u
    double m_pitch_v; // wire pitch v
    double m_pitch_w; // wire pitch w
    double m_ts_width; // time slice width 2 us * 1.6 mm/us ~ 3.2 mm
    double m_angle_u;
    double m_angle_v;
    double m_angle_w;
  public:
    // set defaults
  TPCParams() : m_pitch_u(3)
      , m_pitch_v(3)
      , m_pitch_w(3)
      , m_ts_width(3.2)
      , m_angle_u(1.0472)
      , m_angle_v(-1.0472)
      , m_angle_w(0)
      {};
    
    // set/get u pitches
    void set_pitch_u(double p) { m_pitch_u = p; }
    double get_pitch_u() { return m_pitch_u; }

    // set/get v pitches
    void set_pitch_v(double p) { m_pitch_v = p; }
    double get_pitch_v() { return m_pitch_v; }
    
    // set/get w pitches
    void set_pitch_w(double p) { m_pitch_w = p; }
    double get_pitch_w() { return m_pitch_w; }
    
    double get_pitch(){return (m_pitch_u + m_pitch_v + m_pitch_w)/3.;};
    
    void set_ts_width(double p){ m_ts_width = p;}
    double get_ts_width(){return m_ts_width;}

    void set_angle_u(double p){m_angle_u = p;}
    double get_angle_u(){return m_angle_u;};

    void set_angle_v(double p){m_angle_v = p;}
    double get_angle_v(){return m_angle_v;};

    void set_angle_w(double p){m_angle_w = p;}
    double get_angle_w(){return m_angle_w;};
    

    // etc for other parameters you need
  };
 }
