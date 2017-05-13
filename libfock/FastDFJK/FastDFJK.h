/**
 * Class FastDFJK
 *
 * JK implementation using sieved, threaded
 * density-fitted technology
 */
class FastDFJK : public JK {

protected:

    // => DF-Specific stuff <= //

    /// Auxiliary basis set
    std::shared_ptr<BasisSet> auxiliary_;
    /// PSIO object
    std::shared_ptr<PSIO> psio_;
    /// Cache action for three-index integrals
    std::string df_ints_io_;
    /// Number of threads for DF integrals
    int df_ints_num_threads_;
    /// Condition cutoff in fitting metric, defaults to 1.0E-12
    double condition_;
    /// File number for (Q|mn) tensor
    unsigned int unit_;
    /// Core or disk?
    bool is_core_;
    /// Sieve, must be static throughout the life of the object
    std::shared_ptr<ERISieve> sieve_;
    /// Fitting metric (COULOMB or EWALD) [EWALD is SR]
    std::string metric_;
    /// Ewald metric range parameter
    double theta_;
    /// Geometric atom domain selection algorithm
    std::string domains_;
    /// Flat radius in MHG bump function
    double bump_R0_;
    /// Annihilation radius in MHG bump function
    double bump_R1_;

    // => Required Algorithm-Specific Methods <= //

    /// Do we need to backtransform to C1 under the hood?
    virtual bool C1() const { return true; }
    /// Setup integrals, files, etc
    virtual void preiterations();
    /// Compute J/K for current C/D
    virtual void compute_JK();
    /// Delete integrals, files, etc
    virtual void postiterations();

    /// Common initialization
    void common_init();

    // => Local three-center integrals <= //

    /// Significant atom pairs (reduced triangular indexing A >= B)
    std::vector<std::pair<int, int> > atom_pairs_;
    /// Significant shell pairs (reduced triangular indexing P >= Q)
    std::vector<std::vector<std::pair<int, int> > > shell_pairs_;
    /// Auxiliary basis centers for each atom pair
    std::vector<std::vector<int> > auxiliary_atoms_;
    /// Modified MHG bump function, by atom pair and auxiliary basis center
    std::vector<std::vector<double> > bump_atoms_;
    /// Three-index tensors (pq|A)(A|B)^-1 for each atom pair
    std::vector<std::shared_ptr<Matrix> > Bpq_;

    /// The DF Z operator
    std::shared_ptr<Matrix> Z_;
    /// The DF long-range Z operator
    std::shared_ptr<Matrix> Z_LR_;

    std::shared_ptr<Matrix> build_Z(double omega);
    void build_atom_pairs();
    void build_shell_pairs();
    void build_auxiliary_partition();
    void build_Bpq();
    void bump(std::shared_ptr<Matrix> J,
              const std::vector<double>& bump_atoms,
              const std::vector<int>& auxiliary_atoms,
              bool bump_diagonal);
    void build_J(std::shared_ptr<Matrix> Z,
                 const std::vector<std::shared_ptr<Matrix> >& D,
                 const std::vector<std::shared_ptr<Matrix> >& J);
    void build_K(std::shared_ptr<Matrix> Z,
                 const std::vector<std::shared_ptr<Matrix> >& D,
                 const std::vector<std::shared_ptr<Matrix> >& K);

public:
    // => Constructors < = //

    /**
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     * @param auxiliary auxiliary basis set for this system.
     */
    FastDFJK( std::shared_ptr<BasisSet> primary,
       std::shared_ptr<BasisSet> auxiliary);

    /// Destructor
    virtual ~FastDFJK();

    // => Knobs <= //

    /**
     * Minimum relative eigenvalue to retain in fitting inverse
     * All eigenvectors with \epsilon_i < condition * \epsilon_max
     * will be discarded
     * @param condition minimum relative eigenvalue allowed,
     *        defaults to 1.0E-12
     */
    void set_condition(double condition) { condition_ = condition; }
    /**
     * Metric for FastDF fitting
     * @param metric COULOMB or EWALD,
     *       defaults to COULOMB
     */
    void set_df_metric(const std::string& metric) { metric_ = metric; }
    /**
     * Fitting domain selection algorithm
     * @param domains DIATOMIC, SPHERES
     *       defaults to DIATOMIC
     */
    void set_df_domains(const std::string domains) { domains_ = domains; }
    /**
     * Bump function R0 parameter in a.u. (should be <= R1)
     * @param R0 defaults to 0.0
     */
    void set_df_bump_R0(double R0) { bump_R0_ = R0; }
    /**
     * Bump function R1 parameter in a.u. (should be <= R1)
     * @param R1 defaults to 0.0
     */
    void set_df_bump_R1(double R1) { bump_R1_ = R1; }
    /**
     * Range-Separation parameter for EWALD metric fitting
     * @param theta theta ~ 0 is COULOMB, theta ~ INF is OVERLAP,
     *       defaults to 1.0
     */
    void set_df_theta(double theta) { theta_ = theta; }
    /**
     * Which file number should the (Q|mn) integrals go in
     * @param unit Unit number
     */
    void set_unit(unsigned int unit) { unit_ = unit; }
    /**
     * What action to take for caching three-index integrals
     * @param val One of NONE, LOAD, or SAVE
     */
    void set_df_ints_io(const std::string& val) { df_ints_io_ = val; }
    /**
     * What number of threads to compute integrals on
     * @param val a positive integer
     */
    void set_df_ints_num_threads(int val) { df_ints_num_threads_ = val; }

    // => Accessors <= //

    /**
    * Print header information regarding JK
    * type on output file
    */
    virtual void print_header() const;
};
