/**
 * \file   field_class.hpp
 * \brief  Definition of the class Field
 * \see    field_class.cpp
 * \author Okuto Morikawa
 */
#ifndef FIELD_CLASS_INCLUDED
#define FIELD_CLASS_INCLUDED

#include <iostream>
#include <math.h>
#include <random>
#include "Eigen/LU"

using namespace Eigen;
using std::cout;
using std::endl;
using std::cerr;

/** \brief Imaginary unit */
const std::complex<double> IM(0, 1.0);


/**
 * \enum  Distribution
 * \brief Types of normal distributions
 */
enum class Distribution {
  Zero_conf     /** Zero configuration */,
  Gauss_unit    /** Random device with unit deviation */,
  Gauss_L       /** Random device with deviation \a Li/SQRT2 */,
  Gauss_MT_unit /** Mersenne twistor with unit deviation */,
  Gauss_MT_L    /** Mersenne twistor with deviation \a Li/SQRT2 */
};



/**
 * \class Field
 * \brief Generate normal distributions; Compute convolutions
 * \see   Potential
 * \see   PotentialNR
 * \see   Nicolai
 * \see   Scalar
 * \see   NicolaiSol
 * \see   Distribution
 */
class Field {
protected:
  /** \brief Physical box size, \a N_0=N_1 */
  int Li;
  /** \brief Number of superfields */
  int num_f;
  /**
   * \brief   Superfields
   * \details \a i-th superfield: <em>field(:,i)</em>
   */
  MatrixXcd field;

  /**
   * \fn      void setgauss(const int num_col, const double mean = 0.0, const double dev = 1.0)
   * \details Generate normal distributions 
   *          by random device with unit deviation;
   * \see     Distribution::Gauss_unit
   */
  void setgauss(const int num_col,
				const double mean = 0.0, const double dev = 1.0);
  /**
   * \fn      void setgaussl(const int num_col, const double mean = 0.0)
   * \details Generate normal distributions
   *          by random device with deviation \a Li/SQRT2
   * \see     Distribution::Gauss_L
   */
  void setgaussl(const int num_col,
				 const double mean = 0.0);
  /**
   * \fn      void setgaussmt(const int num_col, const double mean = 0.0, const double dev = 1.0)
   * \details Generate normal distributions
   *          by Mersenne twistor with unit deviation
   * \see     Distribution::Gauss_MT_unit
   */
  void setgaussmt(const int num_col,
				  const double mean = 0.0, const double dev = 1.0);
  /**
   * \fn      void setgaussmtl(const int num_col, const double mean = 0.0)
   * \details Generate normal distributions 
   *          by Mersenne twistor with deviation \a Li/SQRT2
   * \see     Distribution::Gauss_MT_L
   */
  void setgaussmtl(const int num_col,
				   const double mean = 0.0);
  /**
   * \fn      void setconf(const int num_col, const Distribution n, const double mean = 0.0, const double dev = 1.0)
   * \details Generate normal distributions with type Distribution
   * \see     Distribution
   */
  void setconf(const int num_col, const Distribution n,
			   const double mean = 0.0, const double dev = 1.0);
  /**
   * \fn      void setconf(const VectorXcd &v, const int num_field = 1)
   * \details Generate from VectorXcd data
   */
  void setconf(const VectorXcd &v, const int num_field = 1);
  /**
   * \fn      void setconf(const VectorXd &v,  const int num_field)
   * \details Generate from VectorXd data
   */
  void setconf(const VectorXd &v,  const int num_field);
  /** 
   * \fn      void setconf(const MatrixXcd &m)
   * \details Generate from MatrixXcd data
   */
  void setconf(const MatrixXcd &m);
  /**
   * \fn      void setconf(const MatrixXcd &m, const int num_field)
   * \details Generate from MatrixXcd data 
   *          and set number of fields with \a num_field
   */
  void setconf(const MatrixXcd &m, const int num_field);

  /**
   * \fn      MatrixXcd vector2matrix(const int n) const
   * \details From vector-type <em>field(p)</em> 
   *          to matrix-type <em>field(p_0,p_1)</em>
   */
  MatrixXcd vector2matrix(const int n) const;
  /**
   * \fn      MatrixXcd field2matrix(const int n) const
   * \details From vector-type <em>field(p,n)</em> 
   *          to matrix-type <em>field(p-q,n)</em>
   */
  MatrixXcd field2matrix(const int n) const;
  /**
   * \fn      MatrixXcd dfield2matrix(const int n) const
   * \details From vector-type <em>field(p,n)</em> 
   *          to matrix-type <em>field(p,n)*field(q,n)</em>
   */
  MatrixXcd dfield2matrix(const int n) const;

public:
  /** \fn   explicit Field()
   * \brief Constructor of Field */
  explicit Field() {};
  /** 
   * \fn    explicit Field(const int n1, const int n2)
   * \brief Set \a Li=n1 and \a num_f=n2
   * \param n1 Set \a Li
   * \param n2 Set \a num_f
   */
  explicit Field(const int n1, const int n2) {Li = n1; num_f = n2;};
  /**
   * \fn    Field(const int n1, const int n2, const Distribution n3, const double mean = 0.0, const double dev = 1.0)
   * \brief Set \a Li, \a num_f, and \a field
   * \param n1 Set \a Li
   * \param n2 Set \a num_f
   * \param n3 Set normal distribution type
   * \param mean Parameter of normal distribution
   * \param dev  Parameter of normal distribution
   * \see   setconf
   */
  Field(const int n1, const int n2, const Distribution n3,
		const double mean = 0.0, const double dev = 1.0) {
    Li = n1; num_f = n2; setconf(num_f, n3, mean, dev);
  };
  /**
   * \fn    Field(const int n1, const int n2, const int n3, const Distribution n4, const double mean = 0.0, const double dev = 1.0)
   * \brief Set \a Li, \a num_f, and \a field;
   *        \a n3 is 
   * \param n1 Set \a Li
   * \param n2 Set \a num_f
   * \param n3 Set \a <em>field.cols()</em>;
   *        It is possible that <em>n3!=num_f</em>
   * \param n4 Set normal distribution type
   * \param mean Parameter of normal distribution
   * \param dev  Parameter of normal distribution
   * \see   setconf
   */
  Field(const int n1, const int n2, const int n3, const Distribution n4,
		const double mean = 0.0, const double dev = 1.0) {
    Li = n1; num_f = n2; setconf(n3, n4, mean, dev);
  };
  /**
   * \fn    Field(const VectorXcd &v, const int num_field = 1)
   * \brief Set \a Li, \a num_f, and \a field from VectorXcd
   * \see   setconf
   */
  Field(const VectorXcd &v, const int num_field = 1) {setconf(v, num_field);};
  /**
   * \fn    Field(const VectorXd &v, const int num_field)
   * \brief Set \a Li and \a field from VectorXd;
   *        \a num_f!=field.cols()
   * \see   setconf
   */
  Field(const VectorXd &v, const int num_field) {setconf(v, num_field);};
  /**
   * \fn    Field(const MatrixXcd &m)
   * \brief Set \a Li, \a num_f, and \a field from MatrixXcd
   * \see   setconf
   */
  Field(const MatrixXcd &m) {setconf(m);};
  /**
   * \fn      Field(const MatrixXcd &m, const int num_field)
   * \brief   Set \a Li and \a field from MatrixXcd,
   *          but <em>num_f=num_field</em>
   * \details It is possible that <em>num_f!=field.cols()</em>
   * \see     setconf
   */
  Field(const MatrixXcd &m, const int num_field) {setconf(m, num_field);};
  /** \fn virtual ~Field() 
   * \brief Destructor of Field */
  virtual ~Field() {};

  /** \fn static Field Zero(const int n1, const int n2)
   * \see setconf */
  static Field Zero(const int n1, const int n2) {
    return Field {n1, n2, Distribution::Zero_conf}; };
  /** \fn static Field Zero(const int n1, const int n2, const int n3)
   * \see setconf */
  static Field Zero(const int n1, const int n2, const int n3) {
    return Field {n1, n2, n3, Distribution::Zero_conf}; };
  /** \fn static Field Gauss(const int n1, const int n2, const double mean = 0.0, const double dev = 1.0)
   * \see setgauss */
  static Field Gauss(const int n1, const int n2,
					 const double mean = 0.0, const double dev = 1.0) {
    return Field {n1, n2, Distribution::Gauss_unit, mean, dev}; };
  /** \fn static Field Gauss(const int n1, const int n2, const int n3, const double mean = 0.0, const double dev = 1.0) 
   * \see setgauss */
  static Field Gauss(const int n1, const int n2, const int n3,
					 const double mean = 0.0, const double dev = 1.0) {
    return Field {n1, n2, n3, Distribution::Gauss_unit, mean, dev}; };
  /** \fn static Field GaussL(const int n1, const int n2, const double mean = 0.0)
   * \see setgaussl */
  static Field GaussL(const int n1, const int n2,
					  const double mean = 0.0) {
    return Field {n1, n2, Distribution::Gauss_L, mean}; };
  /** \fn static Field GaussL(const int n1, const int n2, const int n3, const double mean = 0.0)
   * \see setgaussl */
  static Field GaussL(const int n1, const int n2, const int n3,
					  const double mean = 0.0) {
    return Field {n1, n2, n3, Distribution::Gauss_L, mean}; };
  /** \fn static Field GaussMT(const int n1, const int n2, const double mean = 0.0, const double dev = 1.0) 
   * \see setgaussmt */
  static Field GaussMT(const int n1, const int n2,
					   const double mean = 0.0, const double dev = 1.0) {
    return Field {n1, n2, Distribution::Gauss_MT_unit, mean, dev}; };
  /** \fn static Field GaussMT(const int n1, const int n2, const int n3, const double mean = 0.0, const double dev = 1.0)
   * \see setgaussmt */
  static Field GaussMT(const int n1, const int n2, const int n3,
					   const double mean = 0.0, const double dev = 1.0) {
    return Field {n1, n2, n3, Distribution::Gauss_MT_unit, mean, dev}; };
  /** \fn static Field GaussMTL(const int n1, const int n2, const double mean = 0.0)
   *\see setgaussmtl */
  static Field GaussMTL(const int n1, const int n2,
						const double mean = 0.0) {
    return Field {n1, n2, Distribution::Gauss_MT_L, mean}; };
  /** \fn static Field GaussMTL(const int n1, const int n2, const int n3, const double mean = 0.0)
   * \see setgaussmtl */
  static Field GaussMTL(const int n1, const int n2, const int n3,
						const double mean = 0.0) {
    return Field {n1, n2, n3, Distribution::Gauss_MT_L, mean}; };
  /** \fn static Field Vector(const VectorXcd &v, const int num_field = 1)
   * \see setconf */
  static Field Vector(const VectorXcd &v, const int num_field = 1) {
    return Field {v, num_field}; };
  /** \fn static Field Matrix(const MatrixXcd &m)
   * \see setconf */
  static Field Matrix(const MatrixXcd &m) {return Field {m};};
  /** \fn static Field Matrix(const MatrixXcd &m, const int num_field) 
   * \see setconf */
  static Field Matrix(const MatrixXcd &m, const int num_field) {
    return Field {m, num_field}; };

  /**
   * \fn      MatrixXcd conf() const
   * \details All superfields
   */
  MatrixXcd conf() const {return field;};
  /**
   * \fn      VectorXcd conf(const int n) const
   * \details \a n-th superfield
   */
  VectorXcd conf(const int n) const {return field.col(n);};
  /**
   * \fn      VectorXd conf_real(const int n) const
   * \details Real part of \a n-th superfield
   */
  VectorXd conf_real(const int n) const {return field.col(n).real();};
  /**
   * \fn      VectorXd conf_imag(const int n) const
   * \details Imaginary part of \a n-th superfield
   */
  VectorXd conf_imag(const int n) const {return field.col(n).imag();};

  /**
   * \fn    virtual bool operator==(const Field &f)
   * \brief Is identical (\a Li, \a num_f) ?
   */
  virtual bool operator==(const Field &f) {return (Li==f.Li) && (num_f==f.num_f);};
  /**
   * \fn    virtual bool operator!=(const Field &f)
   * \brief Not identical (\a Li, \a num_f) ?
   */
  virtual bool operator!=(const Field &f) {return !(*this==f);};

  /**
   * \fn    virtual Field &operator*=(const double n)
   * \brief Multiply by a real number \a n
   */
  virtual Field &operator*=(const double n) { field *= n; return *this; };
  /**
   * \fn    virtual Field &operator/=(const double n)
   * \brief Devide by a real number \a n
   */
  virtual Field &operator/=(const double n) { field *= 1.0/n; return *this;};

  /** 
   * \fn    Field conv(const int n1, const int n2) const
   * \brief Convolution <em>field(:,n1)*field(:,n2)</em>
   * \param n1 \a n1-th superfield
   * \param n2 \a n2-th superfield
   */
  Field conv(const int n1, const int n2) const;
  /**
   * \fn    Field conjconv(const int n) const
   * \brief Convolution <em>field(:,n)*conf(field(:,n))</em>
   * \param n \a n-th superfield
   */
  Field conjconv(const int n) const;
  /**
   * \fn    Field conv_pw(const int pw, const int n = 0) const
   * \brief Convolution <em>field(:,n)*field(:,n)*...*field(:,n)</em>
   * \param pw Convolution times
   * \param n \a n-th superfield
   */
  Field conv_pw(const int pw, const int n = 0) const;

  /**
   * \fn    Field combine_with(const Field &f)
   * \brief Combine with another Field object;
   *        Mutate \a field (\a Li and \a num_f are unchaged)
   *        except for the case that \a field is empty
   */
  Field combine_with(const Field &f);

  /**
   * \fn    virtual void show() const
   * \brief Output \a Li, \a field
   */
  virtual void show() const;
};


#endif
