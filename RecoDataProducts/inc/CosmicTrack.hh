#ifndef RecoDataProducts_CosmicTrack_hh
#define RecoDataProducts_CosmicTrack_hh
//S. Middleton - cosmic track class, will store info of fit result in terms of m and c and errors.




namespace mu2e {

  class CosmicTrack{
  public:

    // Constructors
    CosmicTrack();
    CosmicTrack(double c, double m);
    CosmicTrack(int N, double c, double c_err, double m, double m_err,
               double chisq, double chisq_dof);

    // Destructor
    ~CosmicTrack();

    // Accessors:
   
    int get_N() const { return _Nhits;}
    double get_c() const { return _c; }
    double get_c_err() const { return _c_err; }
    double get_m() const { return _m; }
    double get_m_err() const { return _m_err; }
    double get_chisq() const { return _chisq; }
    double get_chisq_dof() const { return _chisq_dof; }
    //TMatrixD get_coeff_matrix_element(int row, int column) const { return _coeff_matrix[row][column]; }
   
    // Modiffiers:
    void clear();
   
    void set_N(int N){_Nhits = N;}
    void set_c(double c) { _c = c; }
    void set_c_err(double c_err) { _c_err = c_err; }
    void set_m(double m) { _m = m; }
    void set_m_err(double m_err) { _m_err = m_err; }
    void set_chisq(double chisq) { _chisq = chisq; }
    void set_chisq_dof(double chisq_dof) { _chisq_dof = chisq_dof; }
    //void set_coeff_matrix_element(double value, int row, int column) { _coeff_matrix[row][column] = value; }
    
    void set_parameters(int N, double c, double c_err, double m, double m_err, double chisq, double chisq_dof);

  private:
    
    int _Nhits;
    double _c;
    double _c_err;
    double _m;
    double _m_err;
    double _chisq;
    double _chisq_dof; 
   
    
  };
}

#endif

