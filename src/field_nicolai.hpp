/**
 * \file   field_nicolai.hpp
 * \brief  Definition of the classes associated with Nicolai map
 * \see    field_nicolai.cpp
 * \see    field_spt.cpp
 * \author Okuto Morikawa
 */
#ifndef FIELD_NICOLAI_INCLUDED
#define FIELD_NICOLAI_INCLUDED

#include "field_potential.hpp"

#include <iomanip>
#include <fstream>
#include <time.h>   // For test_nr_method //
#include <limits.h> // For test_nr_method //

using std::vector; // For class NicolaiSol //

/** 
 * \brief Identification of scalar solutions
 * \details
 *  If the norm of difference between two Scalars < SOL_ID_MAXVAL,
 *  we regard those configurations as identical ones
 * \see   Scalar
 */
extern const double SOL_ID_MAXVAL;
/**
 * \brief Max error of NR iteration
 * \see   NicolaiSol
 * \see   NicolaiSol::nr_method
 * \see   NicolaiSol::test_nr_method
 */
extern const double MAX_NRERR;
/**
 * \brief   Interruption of NR iteration
 * \details Interrupt NR iteration if nr_err / min(nr_err) >= this parameter
 * \see     NicolaiSol
 * \see     NicolaiSol::nr_method
 * \see     NicolaiSol::test_nr_method
 */
extern const double NR_INTERRUPTION;
/**
 * \brief Standard number of solutions for SuperPotentialType::Custom
 * \see   SuperPotentialType_StdNumSol
 */
extern const int SuperPotentialTypeCustom_StdNumSol;

/**
 * \enum  SuperPotentialType
 * \brief Types of superpotentials
 * \see   field_spt.cpp
 */
enum class SuperPotentialType {
  AlgebraA /** A_n: x^{n+1} */,
  AlgebraD /** D_n: x^{n-1} + x y^2 */,
  AlgebraE /** E_6: x^3 + y^4, E_7: x^3 + x y^3, E_8: x^3 + y^5 */,
  Custom   /** User's setting */
};


/**
 * \class Nicolai
 * \brief Nicolai map; Compute Vector for NR method
 * \see   Field
 * \see   NicolaiSol
 */
class Nicolai:public Field {
protected:
  /** \brief ID number for the Nicolai map */
  int in;

public:
  /** \fn explicit Nicolai()
   * \brief Constructor of Nicolai */
  explicit Nicolai() {};
  /**
   * \fn    explicit Nicolai(const int n1, const int n2, const int n3)
   * \brief Set \a Li, \a num_f, \a in;
   *        Configuration is generated 
   *        by random device with deviation \a Li/SQRT2
   * \param n1 Set \a Li
   * \param n2 Set \a num_f
   * \param n3 Set \a in
   */
  explicit Nicolai(const int n1, const int n2, const int n3) :
    Field(n1, n2, Distribution::Gauss_L) {in = n3;};
  /** \fn virtual ~Nicolai()
   * \brief Destructor of Nicolai */
  virtual ~Nicolai() {};

  /**
   * \fn    virtual bool operator==(const Nicolai &f)
   * \brief Is identical (\a Li, \a num_f) ?
   */
  virtual bool operator==(const Nicolai &f) {return (Li==f.Li) && (num_f == f.num_f);};
  /** 
   * \fn    virtual bool operator!=(const Nicolai &f)
   * \brief Not identical (\a Li, \a num_f) ?
   */
  virtual bool operator!=(const Nicolai &f) {return !(*this==f);};

  /**
   * \fn    VectorXd nr_nvec() const
   * \brief Compute Vector for NR method (Real type)
   * \see   Scalar
   * \see   Scalar::nr_loop
   * \see   PotentialNR
   */
  VectorXd nr_nvec() const;

  /**
   * \fn    VectorXcd nrerr_nvec() const
   * \brief Compute Vector for NR error estimate (complex type)
   * \see   Scalar
   * \see   Scalar::nr_error_vec
   * \see   Scalar::nr_error
   * \see   PotentialNR
   */
  VectorXcd nrerr_nvec() const;

  /**
   * \fn    void nic_output(const VectorXi k, const VectorXd lambda, const int num_nrsol, const VectorXi signs, const int dmp, const double max_err, const SuperPotentialType spt) const
   * \brief Output to file
   * \param k         Power in superpotential
   * \param lambda    Coupling
   * \param num_nrsol Number of "convergent" trials of NR method
   * \param signs     Sign determinant
   * \param dmp       Number of omitting solutions
   *                  because of divergence
   * \param max_err   Maximum error of NR method
   * \param spt       Superpotential type
   */
  void nic_output(const VectorXi k, const VectorXd lambda,
				  const int num_nrsol,
				  const VectorXi signs, const int dmp,
				  const double max_err,
				  const SuperPotentialType spt) const;

  /**
   * \fn    void show() const
   * \brief Output \a Li, \a in, \a field
   */
  void show() const;
};



/**
 * \class Scalar
 * \brief An solution of Nicolai map;
 *        Update to a new solution with NR method;
 *        Compute some types of superpotential; Identify solutions
 * \see   Field
 * \see   Potential
 * \see   PotentialNR
 * \see   NicolaiSol
 * \see   SuperPotentialType
 */
class Scalar:public Field {
  SuperPotentialType spt_type = SuperPotentialType::AlgebraA;
  VectorXi k;
  VectorXd lambda;

  double sol_id_maxval = SOL_ID_MAXVAL;

  /**
   * \fn    Potential superpotential_typealgebraA() const
   * \brief Compute superpotential (class Potential)
   * \see   Potential
   */
  Potential superpotential_typealgebraA() const;
  /**
   * \fn    Potential superpotential_typealgebraD() const
   * \brief Compute superpotential (class Potential)
   * \see   Potential
   */
  Potential superpotential_typealgebraD() const;
  /**
   * \fn    Potential superpotential_typealgebraE() const
   * \brief Compute superpotential (class Potential)
   * \see   Potential
   */
  Potential superpotential_typealgebraE() const;
  /**
   * \fn    Potential superpotential_typecustom() const
   * \brief Compute superpotential (class Potential)
   * \see   Potential
   */
  Potential superpotential_typecustom() const;

  /**
   * \fn    PotentialNR superpotential_nr_typealgebraA() const
   * \brief Compute superpotential (class PotentialNR)
   * \see   PotentialNR
   */
  PotentialNR superpotential_nr_typealgebraA() const;
  /**
   * \fn    PotentialNR superpotential_nr_typealgebraD() const
   * \brief Compute superpotential (class PotentialNR)
   * \see   PotentialNR
   */
  PotentialNR superpotential_nr_typealgebraD() const;
  /**
   * \fn    PotentialNR superpotential_nr_typealgebraE() const
   * \brief Compute superpotential (class PotentialNR)
   * \see   PotentialNR
   */
  PotentialNR superpotential_nr_typealgebraE() const;
  /**
   * \fn    PotentialNR superpotential_nr_typecustom() const
   * \brief Compute superpotential (class PotentialNR)
   * \see   PotentialNR
   */
  PotentialNR superpotential_nr_typecustom() const;

public:
  /** \fn explicit Scalar()
   * \brief Constructor of Scalar */
  explicit Scalar() {};
  /**
   * \fn    explicit Scalar(const int n1, const int n2, const VectorXi n3, const VectorXd n4, const SuperPotentialType t = SuperPotentialType::AlgebraA)
   * \brief Set \a Li, \a num_f, \a k, \a lambda, \a spt_type,
   *        and field configuration generated
   *        by Mersenne twistor with unit deviation
   * \param n1 Set \a Li
   * \param n2 Set \a num_f
   * \param n3 Set \a k
   * \param n4 Set \a lambda
   * \param t  Set \a spt_type
   */
  explicit Scalar(const int n1, const int n2,
				  const VectorXi n3, const VectorXd n4,
				  const SuperPotentialType t
				  = SuperPotentialType::AlgebraA) :
    Field(n1, n2, Distribution::Gauss_MT_unit) {
	k = n3; lambda = n4; spt_type = t;
  };
  /**
   * \fn    explicit Scalar(const VectorXcd &v, const VectorXi n3, const VectorXd n4, const int num_field = 1, const SuperPotentialType t = SuperPotentialType::AlgebraA)
   * \brief Generate from VectorXcd; 
   *        Set \a num_f, \a k, \a lambda, and \a spt_type
   * \param v         Set \a Li and \a field
   * \param n3        Set \a k
   * \param n4        Set \a lambda
   * \param num_field Set \a num_f
   * \param t         Set \a spt_type
   */
  explicit Scalar(const VectorXcd &v,
				  const VectorXi n3, const VectorXd n4,
				  const int num_field = 1,
				  const SuperPotentialType t
				  = SuperPotentialType::AlgebraA) :
    Field(v, num_field) { k = n3; lambda = n4; spt_type = t;};
  /**
   * \fn    explicit Scalar(const Scalar &f, const VectorXcd &v)
   * \brief Set parameters identical to Scalar f;
   *        Configuration is given by VectorXcd
   */
  explicit Scalar(const Scalar &f, const VectorXcd &v) :
    Field(v, f.num_f) {
    spt_type = f.spt_type;
    k = f.k; lambda = f.lambda; spt_type = f.spt_type;
  };
  /**
   * \fn    explicit Scalar(const Scalar &f, const VectorXd &v)
   * \brief Set parameters identical to Scalar f;
   *        Configuration is given by VectorXd
   */
  explicit Scalar(const Scalar &f, const VectorXd &v) :
    Field(v, f.num_f) {
    spt_type = f.spt_type;
    k = f.k; lambda = f.lambda; spt_type = f.spt_type;
  };
  
  /** \fn ~Scalar()
   * \brief Destructor of Scalar */
  ~Scalar() {};

  /**
   * \fn    bool operator==(const Scalar &f)
   * \brief Is identical
   *        (\a Li, \a num_f, \a k, \a lambda, \a spt_type) ?
   */
  bool operator==(const Scalar &f) {
    return (Li==f.Li) && (num_f == f.num_f)
	  && (k==f.k) && (lambda==f.lambda) && (spt_type==f.spt_type);
  };
  /**
   * \fn bool operator!=(const Scalar &f)
   * \brief Not identical
   *        (\a Li, \a num_f, \a k, \a lambda, \a spt_type) ?
   */
  bool operator!=(const Scalar &f) { return !(*this==f); };

  /**
   * \fn      bool is_identical(const Scalar &f) const
   * \brief   Identify two solutions
   * \param   f Another configuration to be compared with this
   * \return  if (identical) true;
   *          else false
   * \details Threshold: <em>sol_id_maxval</em>
   */
  bool is_identical(const Scalar &f) const;

  /**
   * \fn    VectorXcd nrerr_svec() const
   * \brief Compute Vector for NR error (complex type)
   * \see   NR_error
   * \see   Nicolai
   * \see   PotentialNR
   */
  VectorXcd nrerr_svec() const;

  /**
   * \fn    Potential superpotential() const
   * \brief Compute superpotential (class Potential)
   * \see   Potential
   */
  Potential superpotential() const;
  /**
   * \fn    PotentialNR superpotential_nr() const
   * \brief Compute superpotential (class PotentialNR)
   * \see   PotentialNR
   */
  PotentialNR superpotential_nr() const;

  /**
   * \fn     Scalar nr_loop (const Nicolai &nic)
   * \brief  An iteration of NR method;
   *         Compute LU decompositon; Overwrite own members
   * \param  nic Nicolai map
   * \return New configuration of scalar solution
   */
  Scalar nr_loop (const Nicolai &nic);

  /**
   * \fn     VectorXcd nr_error_vec (const Nicolai &nic) const
   * \brief  Compute error of NR method for each momentum
   * \param  nic Nicolai map
   * \return Error estimation of NR method for each momentum
   */
  VectorXcd nr_error_vec (const Nicolai &nic) const;
  /**
   * \fn     double nr_error (const Nicolai &nic) const
   * \brief  Compute error of NR method
   * \param  nic Nicolai map
   * \return Error estimation of NR method; nr_error_vec.norm()
   */
  double nr_error (const Nicolai &nic) const {
    VectorXcd er = nr_error_vec(nic);
    return er.norm();
  };

  /**
   * \fn    void scl_output(const int in, const int num_sol, const int i, const int sign, const double err)
   * \param in      ID number of the Nicolai map
   * \param num_sol Number of solutions
   * \param i       ID number of the solution
   * \param sign    Sign determinant
   * \param err     Error of NR method
   * \brief Output to file
   */
  void scl_output(const int in, const int num_sol, const int i,
				  const int sign, const double err) const;

  /**
   * \fn    void show() const
   * \brief Output \a Li, \a k, \a lambda, \a field
   */
  void show() const;
};



/**
 * \class NicolaiSol
 * \brief Execute the Newton--Raphson method;
 *        Combine solutions; Obtain sign det for each Scalar
 * \see   Scalar
 */
class NicolaiSol:public Nicolai {
  vector<Scalar> fs;
  vector<double> errs;

  double nr_interruption = NR_INTERRUPTION;

public:
  /** \fn explicit NicolaiSol()
   * \brief Constructor of NicolaiSol */
  explicit NicolaiSol() {};
  /**
   * \fn    explicit NicolaiSol(const int n1, const int n2, const int n3)
   * \brief Set members of the class Nicolai:
   * \param n1 Set \a Li
   * \param n2 Set \a num_f
   * \param n3 Set in
   */
  explicit NicolaiSol(const int n1, const int n2, const int n3) :
	Nicolai(n1, n2, n3) {};
  /** \fn ~NicolaiSol ()
   * \brief Destructor of NicolaiSol */
  ~NicolaiSol () {};

  /**
   * \fn      bool operator==(const NicolaiSol &f)
   * \brief   Is identical (\a Li, \a num_f, field) ?
   * \details Note that the operator compares <em>field</em>s
   */
  bool operator==(const NicolaiSol &f) {
	MatrixXcd err_vec = (field-f.field) / field.norm();
	double err = err_vec.norm();
	if (err < 1.0e-15)
	  return (Li==f.Li) && (num_f == f.num_f);
	else return false;
  };
  /** 
   * \fn      bool operator!=(const NicolaiSol &f)
   * \brief   Not identical (\a Li, \a num_f, field) ?
   * \details Note that the operator compares <em>field</em>s
   */
  bool operator!=(const NicolaiSol &f) {return !(*this==f);};

  /**
   * \fn    NicolaiSol &add_sol(const Scalar &f, const double &err)
   * \brief Add a new solution; Ignore identical one
   */
  NicolaiSol &add_sol(const Scalar &f, const double &err);
  /**
   * \fn    NicolaiSol &add_sol(const NicolaiSol &sol)
   * \brief Add a new solution; Ignore identical ones
   */
  NicolaiSol &add_sol(const NicolaiSol &sol);

  /** \fn   int num_sol() const
   * \brief Number of solutions */
  int num_sol() const {return fs.size();};
  /** \fn   double max_err() const
   * \brief Maximum error of NR method */
  double max_err() const {return *max_element(errs.begin(), errs.end());};

  /**
   * \fn     double nr_method (const VectorXi k, const VectorXd lambda, const int num_nrsol, const int LOOP)
   * \brief  Newton--Raphson method
   * \param  k         Power in superpotential
   * \param  lambda    Coupling
   * \param  num_nrsol Number of "convergent" trials of NR method
   * \param  LOOP      Maximum number of iteration
   * \param  t         Type of SuperPotential
   * \return Number of omitting initial configuraiotns
   *         because of divergence
   */
  double nr_method (const VectorXi k, const VectorXd lambda,
					const int num_nrsol, const int LOOP,
					const SuperPotentialType t
					= SuperPotentialType::AlgebraA);

  /**
   * \fn    void test_nr_method (const VectorXi k, const VectorXd lambda, const int LOOP)
   * \brief Observe the numerical convergence 
   *        of NR method for each system
   * \param k      Power in superpotential
   * \param lambda Coupling
   * \param LOOP   Maximum number of iteration
   * \param TRIAL  Maximum number of trials
   * \param t      Type of SuperPotential
   */
  void test_nr_method (const VectorXi k, const VectorXd lambda,
					   const int LOOP, const int TRIAL,
					   const SuperPotentialType t
					   = SuperPotentialType::AlgebraA);

  /**
   * \fn    int sign_det(const int n) const
   * \brief Sign determinant (Jacobian)
   */
  int sign_det(const int n) const {
    return fs[n].superpotential().sign_det();
  };

  /**
   * \fn    void phi_output(const int in, const VectorXi signs) const
   * \brief Output to file
   */
  void phi_output(const int in, const VectorXi signs) const;

  /**
   * \fn    void show() const
   * \brief Output Scalars
   */
  void show() const;
};

/**
 * \fn    int SuperPotentialType_StdNumSol (const VectorXi k, const SuperPotentialType t)
 * \brief Standard number of solutions for each SuperPotentialType
 */
extern int SuperPotentialType_StdNumSol (const VectorXi k, const SuperPotentialType t);

#endif
