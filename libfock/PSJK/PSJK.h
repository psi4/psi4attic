
/**
 * Class PSJK
 *
 * JK implementation using sieved, threaded,
 * range-separated PS technology
 */
class PSJK : public JK {

protected:

    /// Options reference (needed to build grid)
    Options& options_;
    /// Dealiasing basis (if needed)
    std::shared_ptr<BasisSet> dealias_;
    /// Number of threads for three-center integrals
    int df_ints_num_threads_;
    /// File number for (Q|mn) tensor
    unsigned int unit_;
    /// Range separation for integrand smoothing
    double theta_;
    /// QUADRATURE, RENORMALIZATION, or DEALIASING
    std::string dealiasing_;
    /// PSIO object
    std::shared_ptr<PSIO> psio_;
    /// Sieve, must be static throughout the life of the object
    std::shared_ptr<ERISieve> sieve_;
    /// Q_m^P matrix
    SharedMatrix Q_;
    /// R_m^P matrix
    SharedMatrix R_;
    /// Grid definition [P x [x y z w]]
    SharedMatrix grid_;
    /// 4-center integrators
    std::vector<std::shared_ptr<TwoBodyAOInt> > ints_4c_;
    /// Q R D (for J)
    SharedVector d_;
    /// R D (for K)
    SharedMatrix V_;
    /// V A (for K)
    SharedMatrix W_;
    /// Temporary triangular J
    SharedVector J_temp_;

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

    // => Magic <= //
    void build_QR();
    void build_Amn_disk(double theta, const std::string& entry);
    void block_J(double** Qmnp, int Pstart, int nP, const std::vector<SharedMatrix>& J);
    void block_K(double** Qmnp, int Pstart, int nP, const std::vector<SharedMatrix>& K);
    void build_JK_SR();
    void build_JK_LR();

    void build_JK_debug(const std::string& op = "", double theta = 0.0);

    int max_rows();

public:

    // => Constructors < = //

    /**
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     * @param options, Options reference used to build grid
     */
    PSJK(std::shared_ptr<BasisSet> primary,
        Options& options);

    /// Destructor
    virtual ~PSJK();

    // => Knobs <= //

    /**
     * What number of threads to compute integrals on?
     * @param val a positive integer
     */
    void set_df_ints_num_threads(int val) { df_ints_num_threads_ = val; }
    /**
     * What value of range-separation parameter for integral smoothing?
     * @param theta a positive double
     */
    void set_theta(double theta) { theta_ = theta; }
    /**
     * How to handle the renormalization or dealiasing?
     * @param type QUADRATURE, RENORMALIZATION, or DEALIASING
     */
    void set_dealiasing(const std::string& dealiasing) { dealiasing_ = dealiasing; }
    /**
     * Custom dealias basis
     * @param dealias, new dealias basis
     */
    void set_dealias_basis(std::shared_ptr<BasisSet> dealias) { dealias_ = dealias; }
    /**
     * Which file number should the (Q|mn) integrals go in
     * @param unit Unit number
     */
    void set_unit(unsigned int unit) { unit_ = unit; }

    // => Accessors <= //

    /**
    * Print header information regarding JK
    * type on output file
    */
    virtual void print_header() const;

};
