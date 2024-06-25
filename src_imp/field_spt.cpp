/**
 * \file   field_spt.cpp
 * \brief  Definition of the methods, Scalar::superpotential and Scalar::superpotential_nr
 *         in the case of SuperPotentialType::Custom
 * \see    field_nicolai.hpp
 * \author Okuto Morikawa
 */
#include "field_nicolai.hpp"

/**
 * \brief Standard number of solutions for SuperPotentialType::Custom
 * \see   SuperPotentialType_StdNumSol
 */
const int SuperPotentialTypeCustom_StdNumSol = 8;


Potential Scalar::superpotential_typecustom () const {
  //
  // lambda0 X^3 + lambda1 Y^3 + lambda2 Z^3 + lambda3 X Y Z
  //
  try {
	if (num_f != 3)
	  throw "ERROR<Scalar::superpotential> Case SuperPotentialType::Custom: number of fields should be 3.";
	if (lambda.size() != 4)
	  throw "ERROR<Scalar::superpotential> Case SuperPotentialType::Custom: size of VectorXd for couplings should be 4.";
	//
	Field X1 = conv_pw(1, 0); X1 *= 2.0 * lambda(0);
	Field Y1 = conv_pw(1, 1); Y1 *= 2.0 * lambda(1);
	Field Z1 = conv_pw(1, 2); Z1 *= 2.0 * lambda(2);

	Field Xlm = conv_pw(1, 0); Xlm *= lambda(3);
	Field Ylm = conv_pw(1, 1); Ylm *= lambda(3);
	Field Zlm = conv_pw(1, 2); Zlm *= lambda(3);

	Field res = X1;
	res.combine_with(Zlm); res.combine_with(Ylm);
	res.combine_with(Y1);  res.combine_with(Xlm);
	res.combine_with(Z1);
	return Potential {res.conf(), num_f};
	//
  } catch (const char *str) {
	cerr << str << endl; return Potential {Li, num_f};
  }
}

PotentialNR Scalar::superpotential_nr_typecustom () const {
  //
  // lambda0 X^3 + lambda1 Y^3 + lambda2 Z^3 + lambda3 X Y Z
  //
  try {
	if (num_f != 3)
	  throw "ERROR<Scalar::superpotential_nr> Case SuperPotentialType::Custom: number of fields should be 3.";
	if (lambda.size() != 4)
	  throw "ERROR<Scalar::superpotential_nr> Case SuperPotentialType::Custom: size of VectorXd for couplings should be 4.";
	//
	Field X1 = conv_pw(1, 0); X1 *= 2.0 * lambda(0);
	Field Y1 = conv_pw(1, 1); Y1 *= 2.0 * lambda(1);
	Field Z1 = conv_pw(1, 2); Z1 *= 2.0 * lambda(2);
	Field Xlm = conv_pw(1, 0); Xlm *= lambda(3);
	Field Ylm = conv_pw(1, 1); Ylm *= lambda(3);
	Field Zlm = conv_pw(1, 2); Zlm *= lambda(3);
	Field res1 = X1;
	res1.combine_with(Zlm); res1.combine_with(Ylm);
	res1.combine_with(Y1);  res1.combine_with(Xlm);
	res1.combine_with(Z1);

	Field X2yz;
	{
	  Field X2 = conv_pw(2, 0); X2 *= lambda(0);
	  Field YZ = conv(1, 2); YZ *= lambda(3);
	  MatrixXcd m = X2.conf() + YZ.conf();
	  X2yz = Field {m};
	}
	Field Y2xz;
	{
	  Field Y2 = conv_pw(2, 1); Y2 *= lambda(1);
	  Field XZ = conv(0, 2); XZ *= lambda(3);
	  MatrixXcd m = Y2.conf() + XZ.conf();
	  Y2xz = Field {m};
	}
	Field Z2xy;
	{
	  Field Z2 = conv_pw(2, 2); Z2 *= lambda(2);
	  Field XY = conv(0, 1); XY *= lambda(3);
	  MatrixXcd m = Z2.conf() + XY.conf();
	  Z2xy = Field {m};
	}
	Field res2 = X2yz;
	res2.combine_with(Y2xz); res2.combine_with(Z2xy);
	MatrixXcd m_res2 = res2.conf();
	return PotentialNR {num_f, res1.conf(), m_res2, 2.0 * m_res2};
	//
  } catch (const char *str) {
	cerr << str << endl; return PotentialNR {Li, num_f};
  }
}
