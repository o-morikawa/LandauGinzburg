/**
 * \file   field_potential.hpp
 * \brief  Definition of the class Potential
 * \see    field_potential.cpp
 * \author Okuto Morikawa
 */
#ifndef FIELD_POTENTIAL_INCLUDED
#define FIELD_POTENTIAL_INCLUDED

#include "field_class.hpp"

/** 
 * \namespace MatrixComp
 * \brief     Mainly functions to compute Jacobian for uni-superfield
 * \see       Potential
 * \see       PotentialNR
 */
namespace MatrixComp {
  /**
   * \fn    MatrixXcd ResizeMatrix (const MatrixXcd &m, const int rw, const int cl)
   * \brief Output m as MatrixXcd(rw, cl) with zeros
   * \see   PotentialNR::PotentialNR
   */
  extern MatrixXcd ResizeMatrix (const MatrixXcd &m,
								 const int rw, const int cl);
  /**
   * \fn    MatrixXcd DeleteZeroLine (const MatrixXcd &m)
   * \brief Delete a row and a column with \a p=0, \a q=0
   * \see   Potential::unifield_tilde_type
   * */
  extern MatrixXcd DeleteZeroLine (const MatrixXcd &m);
  /** 
   * \fn    double DiffLU (const MatrixXcd &m, const PartialPivLU<MatrixXcd> &lu)
   * \brief Difference between a matrix and its LU decomposition 
   * \see   MatrixComp::SignDet
   */
  extern double DiffLU (const MatrixXcd &m,
						const PartialPivLU<MatrixXcd> &lu);
  /**
   * \fn    int SignDet (const MatrixXcd &m, const int debug = 0)
   * \brief Sign of determinant with LU Decomposition
   * \see   MatrixComp::DiffLU
   */
  extern int SignDet (const MatrixXcd &m, const int debug = 0);
};



/**
 * \class   Potential
 * \brief   Compute Jacobian and its sign determinant
 * \details \a field:
 *          <em>partial_i partial_j W({A}) (p) with i <= j</em>;
 *          Then, <em>num_f!=field.cols()</em>,
 *          but <em>n2(n2+1)/2==field.cols()</em>
 * \see     Field
 * \see     PotentialNR
 * \see     MatrixComp
 */
class Potential :public Field {
// protected:
  /**
   * \fn      void setpotential(const MatrixXcd &m)
   * \brief   Set \a Li and \a field from MatrixXcd,
   *          but \a num_f such as
   *          <em>num_f(num_f+1)/2==field.cols()</em>
   * \see     setconf
   */
  void setpotential(const MatrixXcd &m);
  /**
   * \fn      void setpotential(const MatrixXcd &m, const int num_field)
   * \brief   Set \a Li and \a field from MatrixXcd,
   *          but <em>num_f=num_field</em>
   * \details <em>num_f(num_f+1)/2==field.cols()</em>
   * \see     setconf
   */
  void setpotential(const MatrixXcd &m, const int num_field);

  /** 
   * \fn    MatrixXcd jacobian() const
   * \brief Compute Jacobian
   * \see   sign_det
   */
  MatrixXcd jacobian() const;

  /** 
   * \fn    MatrixXcd unifield_tilde_type() const
   * \brief Computation of Jacobian for uni-superfield
   * \see   MatrixComp::DeleteZeroLine
   */
  MatrixXcd unifield_tilde_type() const;
  /**
   * \fn    MatrixXcd unifield_jacobian() const
   * \brief Computation of Jacobian for uni-superfield
   */
  MatrixXcd unifield_jacobian() const;
  /**
   * \fn    MatrixXcd dualfield_jacobian() const
   * \brief Computation of Jacobian for dual-superfield
   */
  MatrixXcd dualfield_jacobian() const;
  /**
   * \fn    MatrixXcd triplefield_jacobian1() const
   * \brief Computation of Jacobian for triple-superfield
   */
  MatrixXcd triplefield_jacobian1() const;
  /**
   * \fn    MatrixXcd triplefield_jacobian2() const
   * \brief Computation of Jacobian for triple-superfield
   */
  MatrixXcd triplefield_jacobian2() const;
  /**
   * \fn    MatrixXcd multifield_jacobian() const
   * \brief Computation of Jacobian for multi-superfield
   */
  MatrixXcd multifield_jacobian() const;

public:
  /** \fn explicit Potential()
   * \brief   Constructor of Potential */
  explicit Potential() {};
  /**
   * \fn    explicit Potential(const int n1, const int n2, const Distribution n3 = Distribution::Zero_conf)
   * \brief Set \a Li, \a num_f, and \a field
   * \param n1 Set \a Li
   * \param n2 Set \a num_f; <em>n2(n2+1)/2==field.cols()</em>
   * \param n3 Set normal distribution type
   * \see   setconf
   */
  explicit Potential(const int n1, const int n2,
					 const Distribution n3 = Distribution::Zero_conf)
    : Field(n1, n2) {
    int num_p = int(num_f*(num_f+1)/2); setconf(num_p, n3); };
  /**
   * \fn    explicit Potential(const MatrixXcd &m)
   * \brief Set \a Li and \a field from MatrixXcd,
   *        but \a num_f such as
   *        <em>num_f(num_f+1)/2==field.cols()</em>
   * \see   setpotential
   */
  explicit Potential(const MatrixXcd &m) {
	setpotential(m);
  };
  /**
   * \fn      explicit Potential(const MatrixXcd &m, const int num_field)
   * \brief   Set \a Li and \a field from MatrixXcd,
   *          but <em>num_f=num_field</em>
   * \details <em>num_f(num_f+1)/2==field.cols()</em>
   * \see     setpotential
   */
  explicit Potential(const MatrixXcd &m, const int num_field) {
	setpotential(m, num_field);
  };
  /** \fn virtual ~Potential()
   * \brief Destructor of Potential */
  virtual ~Potential() {};

  /** 
   * \fn    virtual bool operator==(const Potential &p) 
   * \brief Is identical (\a Li, \a num_f) ?
   */
  virtual bool operator==(const Potential &p) {return (Li==p.Li) && (num_f == p.num_f);};
  /**
   * \fn    virtual bool operator!=(const Potential &p)
   * \brief Not identical (\a Li, \a num_f) ?
   */
  virtual bool operator!=(const Potential &p) {return !(*this==p);};

  /** 
   * \fn    int sign_det(const int debug = 0)
   * \brief Sign determinant with a fast algorithm (<em>num_f < 4</em>)
   * \details
   *        Probability of faster cases than Potential::sign_det4multifield:
   *        <em>num_f=1</em>: 99.9\%,
   *        <em>num_f=2</em>: 99\%
   *        <em>num_f=3</em>: 97-98\%
   * \param debug Debugging mode on(0)/off(1)
   * \see   MatrixComp::SignDet
   */
  int sign_det(const int debug = 0) const;
  /** 
   * \fn    int sign_det4multifield(const int debug = 0)
   * \brief Sign determinant for mutli-superfield
   * \param debug Debugging mode on(0)/off(1)
   * \see   MatrixComp::SignDet
   */
  int sign_det4multifield (const int debug = 0) const;
  /**
   * \fn    virtual void show() const
   * \brief Output \a Li, \a field
   */
  virtual void show() const;
};



/**
 * \class   PotentialNR
 * \brief   Compute Matrix and Vector for NR method
 * \details \a field: 
 *          <em>partial_i partial_j W({A}) (p) with i <= j</em>
 *
 *          \a field2:
 *          <em>partial_i W({A}) (p)</em>
 *
 *          \a field3:
 *          <em>sum_j A_j partial_i partial_j W({A}) (p)</em>
 * \see     Field
 * \see     Potential
 * \see     Nicolai
 * \see     Scalar
 * \see     Scalar::nr_loop
 */
class PotentialNR: public Potential {
  /** \details <em>partial_i W({A}) (p)</em> */
  MatrixXcd field2;
  /** \details <em>sum_j A_j partial_i partial_j W({A}) (p)</em> */
  MatrixXcd field3;

public:
  /** \fn explicit PotentialNR()
   * \brief   Constructor of PotentialNR */
  explicit PotentialNR() {};
  /** 
   * \fn    explicit PotentialNR(const int n1, const int n2, const Distribution n3 = Distribution::Zero_conf)
   * \brief Set \a Li, \a num_f, and \a field-field2-field3 with zeros
   * \param n1 Set \a Li
   * \param n2 Set \a num_f;
   *        <em>n2(n2+1)/2==field.cols()</em>;
   *        <em>n2==field2.cols()==field3.cols()</em>
   * \param n3 Set normal distribution type
   */
  explicit PotentialNR(const int n1, const int n2,
					   const Distribution n3 = Distribution::Zero_conf)
    : Potential(n1, n2, n3) {
    int N = (Li+1)*(Li+1);
    field2 = MatrixXcd::Zero(N, num_f);
    field3 = MatrixXcd::Zero(N, num_f);
  };
  /** 
   * \fn    explicit PotentialNR(const int num_field, const MatrixXcd &m1, const MatrixXcd &m2, const MatrixXcd &m3)
   * \brief Set \a Li, \a num_f, and \a field-field2-field3
   * \param num_field Set \a num_f;
   * \param m1        Set \a field;
   * \param m2        Set \a field2;
   * \param m3        Set \a field3;
   */
  explicit PotentialNR(const int num_field,
					   const MatrixXcd &m1,
					   const MatrixXcd &m2, const MatrixXcd &m3)
    : Potential(m1, num_field) {
	int N = m1.rows();
	try {
	  if (m2.rows()!=N || m3.rows()!=N
		  || m2.cols()!=num_f || m3.cols()!=num_f) throw 1;
	  field2 = m2; field3 = m3;
	} catch (const int i) {
	  field2 = MatrixComp::ResizeMatrix(m2, N, num_f);
	  field3 = MatrixComp::ResizeMatrix(m3, N, num_f);
	}
  };
  /** \fn ~PotentialNR()
   * \brief Destructor of PotentialNR */
  ~PotentialNR() {};

  /**
   * \fn    virtual bool operator==(const PotentialNR &p) 
   * \brief Is identical (\a Li, \a num_f) ?
   */
  bool operator==(const PotentialNR &p) {return (Li==p.Li) && (num_f == p.num_f);};
  /**
   * \fn    virtual bool operator!=(const PotentialNR &p)
   * \brief Not identical (\a Li, \a num_f) ? 
   */
  bool operator!=(const PotentialNR &p) {return !(*this==p);};

  /**
   * \fn    VectorXd nr_pvec() const
   * \brief Compute Vector for NR method (Real type)
   * \see   Nicolai
   * \see   Scalar
   * \see   Scalar::nr_loop
   */
  VectorXd nr_pvec() const;
  /**
   * \fn    MatrixXd nr_mat() const
   * \brief Compute Matrix for NR method (Real type)
   * \see   Nicolai
   * \see   Scalar
   * \see   Scalar::nr_loop
   */
  MatrixXd nr_mat() const;

  /**
   * \fn    VectorXcd nrerr_pvec() const
   * \brief Compute Vector for NR error estimate (complex type)
   * \see   Nicolai
   * \see   Scalar
   * \see   Scalar::nr_error
   */
  VectorXcd nrerr_pvec() const;

  /**
   * \fn    void show() const
   * \brief Output \a Li, \a field, \a field2, \a field3
   */
  void show() const;
};

#endif
