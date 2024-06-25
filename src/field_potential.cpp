/**
 * \file   field_potential.cpp
 * \brief  Definition of the class Potential
 * \see    field_potential.hpp
 * \author Okuto Morikawa
 */
#include "field_potential.hpp"


MatrixXcd MatrixComp::ResizeMatrix (const MatrixXcd &m,
									const int rw, const int cl) {
  MatrixXcd res = MatrixXcd::Zero(rw, cl);
  int r2 = m.rows(); int c2 = m.cols();
  if (r2 >= rw) {
	if (c2 >= cl)
	  res                  = m.block(0,0,rw,cl);
	else
	  res.block(0,0,rw,c2) = m.block(0,0,rw,c2);
  } else {
	if (c2 >= cl)
	  res.block(0,0,r2,cl) = m.block(0,0,r2,cl);
	else
	  res.block(0,0,r2,c2) = m;
  }
  return res;
}

MatrixXcd MatrixComp::DeleteZeroLine (const MatrixXcd &m) {
  int r = m.rows();
  int c = m.cols();
  try {
    if (r % 2 == 0 || c % 2 == 0)
      throw "ERROR<DeleteZeroLine> Size of rows() or cols() is not ODD.";
    int Qr = r / 2;
    int Qc = c / 2;
    MatrixXcd T(r - 1, c - 1);
    T.block(0,  0,  Qr, Qc) = m.block(0,    0,    Qr, Qc);
    T.block(0,  Qc, Qr, Qc) = m.block(0,    Qc+1, Qr, Qc);
    T.block(Qr, 0,  Qr, Qc) = m.block(Qr+1, 0,    Qr, Qc);
    T.block(Qr, Qc, Qr, Qc) = m.block(Qr+1, Qc+1, Qr, Qc);
    return T;
  } catch (const char *str) {cerr << str << endl; return MatrixXcd::Zero(r-1, c-1);}
}
double MatrixComp::DiffLU (const MatrixXcd &m,
						   const PartialPivLU<MatrixXcd> &lu) {
  int N = m.rows();
  try {
    if (N != m.cols())
      throw "ERROR<DiffLU> Not a square matrix.";
    MatrixXcd P, L, U;
    P = lu.permutationP();
    L = lu.matrixLU().triangularView<StrictlyLower>();
    L += MatrixXcd::Identity(N, N);
    U = lu.matrixLU().triangularView<Upper>();
    MatrixXcd Prod; Prod = P.inverse() * L * U;
    MatrixXcd Diff; Diff = m - Prod;
    return Diff.norm();
  } catch (const char *str) {cerr << str << endl; return -1.0;}
}
int MatrixComp::SignDet (const MatrixXcd &m, const int debug) {
  try {
    if (m.rows() != m.cols())
      throw "ERROR<SignDet> Not a square matrix.";
    PartialPivLU<MatrixXcd> lu(m);
    VectorXcd diag; diag = lu.matrixLU().diagonal();
    std::complex<double> phase = 1;
    for (int i = 0; i < diag.size(); ++i) {
      phase *= diag(i) / std::abs(diag(i));
    }
    if ( lu.permutationP().determinant() < 0 ) phase *= -1;
    int sign = 1;
    if (phase.real() < 0) sign *= -1;
    if (debug != 0) {
      double diff = MatrixComp::DiffLU(m, lu);
      VectorXd im = lu.matrixLU().diagonal().imag();
      double im_norm = im.norm();
      double im_max = 0;
      for (int i = 0; i < im.size(); ++i) {
		if (im_max < std::abs(im(i))) im_max = std::abs(im(i));
      }
      cerr << "PartialPiv: " << "Difference   " << diff    << endl
		   << "            " << "Imag  Norm   " << im_norm << endl
		   << "            " << "      Max    " << im_max  << endl
		   << "            " << "Phase of Det " << phase   << endl;
    }
    return sign;
  } catch (const char *str) {cerr << str << endl; return 0;}
}



void Potential::setpotential(const MatrixXcd &m) {
  int N = m.rows(); int mcol = m.cols();
  int num_field; int num_col;
  for (int i = 1; i <= mcol; ++i) {
	int tmp = int(i * (i + 1) / 2);
	if (tmp >= mcol) {num_field = i; num_col = tmp; break;}
  }
  if (num_col == mcol) setconf(m, num_field);
  else {
	MatrixXcd mbig = MatrixXcd::Zero(N, num_col);
	mbig.block(0,0, N, mcol) = m;
	setconf(mbig, num_field);
  }
  return;
}
void Potential::setpotential(const MatrixXcd &m, const int num_field) {
  if (num_field < 1) {setpotential(m); return;}
  int num_col = int(num_field*(num_field+1)/2);
  int mcol = m.cols();
  try {
	if (num_col > mcol)      throw 1;
	else if (num_col < mcol) throw 2;
	setconf(m, num_field);
	return;
  } catch (const int i) {
	int N = m.rows();
	int num_field2; int num_col2;
	switch (i) {
	case 1:
	  num_field2 = num_field; num_col2 = num_col;
	  break;
	case 2:
	  for (int j = num_field; j <= mcol; ++j) {
		int tmp = int(j * (j+1) / 2);
		if (tmp >= mcol) {num_field2 = j; num_col2 = tmp; break;}
	  }
	  break;
	}
	MatrixXcd mbig = MatrixXcd::Zero(N, num_col2);
	mbig.block(0,0, N, mcol) = m;
	setconf(mbig, num_field2);
	return;
  }
}

MatrixXcd Potential::unifield_tilde_type () const {
  int N = (Li+1) * (Li+1);
  try {
    if (N % 2 == 0)
      throw "ERROR<unifield_tilde_type> Size of Vector is not ODD.";
    if (num_f != 1)
      throw "ERROR<unifield_tilde_type> Number of superfields is not unity.";
    MatrixXcd wm; wm = field2matrix(0) - dfield2matrix(0);
    return MatrixComp::DeleteZeroLine(wm);
  } catch (const char *str) {cerr << str << endl; return MatrixXcd::Zero(N-1, N-1);}
}

MatrixXcd Potential::unifield_jacobian () const {
  int N = (Li+1) * (Li+1);
  try {
    if (N % 2 == 0)
      throw "ERROR<unifield_jacobian> Size of Vector is not ODD.";
    if (num_f != 1)
      throw "ERROR<unifield_jacobian> Number of superfields is not unity.";
    int L = Li+1;
    MatrixXcd wt = unifield_tilde_type();
    MatrixXcd p, pi;
    p = MatrixXcd::Zero(N-1, N-1);
    pi = MatrixXcd::Zero(N-1, N-1);
    for (int i = 0; i < N/2; ++i) {
      p(i, i) = 2.0 * M_PI * double(L-1)
		* (double(L/2 - i/L) + IM * double(L/2 - i%L));
      p(N-2 - i, N-2 - i) = - p(i, i);
      pi(i, i) = p(i, i) / std::norm(p(i,i));
      pi(N-2 - i, N-2 - i) = - pi(i,i);
    }
    return p + wt * pi * wt.adjoint();
  } catch (const char *str) {cerr << str << endl; return MatrixXcd::Zero(N-1, N-1);}
}
MatrixXcd Potential::dualfield_jacobian () const {
  int N = (Li+1) * (Li+1);
  try {
    if (N % 2 == 0)
      throw "ERROR<dualfield_jacobian> Size of Vector is not ODD.";
    if (num_f != 2)
      throw "ERROR<dualfield_jacobian> Number of superfields is not 2.";
	int L = Li+1;
	MatrixXcd w12; w12 = field2matrix(1);
	MatrixXcd m1 = MatrixXcd::Zero(2*N, 2*N);
	m1.block(0,0, N,N) = w12; m1.block(N,N, N,N) = w12.adjoint();
	MatrixXcd m2 = MatrixXcd::Zero(2*N, 2*N);
	{
	  MatrixXcd p = MatrixXcd::Zero(N, N);
	  for (int i = 0; i < N; ++i) {
		p(i, i) = 2.0 * M_PI * double(Li) * (IM * double(L/2 - i/L) + double(L/2 - i%L));
	  }
	  MatrixXcd m21 = MatrixXcd::Zero(2*N, 2*N);
	  {
		MatrixXcd w22; w22 = field2matrix(2);
		m21.block(0,0, N,N) = w22; m21.block(N,N, N,N) = w22.adjoint();
		m21.block(N,0, N,N) = p;   m21.block(0,N, N,N) = - p.adjoint();
	  }
	  MatrixXcd m22 = MatrixXcd::Zero(2*N, 2*N);
	  {
		PartialPivLU<MatrixXcd> lu(w12);
		MatrixXcd m_inv; m_inv = lu.inverse();
		m22.block(0,0, N,N) = m_inv; m22.block(N,N, N,N) = m_inv.adjoint();
	  }
	  MatrixXcd m23 = MatrixXcd::Zero(2*N, 2*N);
	  {
		MatrixXcd w11; w11 = field2matrix(0);
		m23.block(0,0, N,N) = w11; m23.block(N,N, N,N) = w11.adjoint();
		m23.block(N,0, N,N) = p;   m23.block(0,N, N,N) = - p.adjoint();
	  }
	  m2 = m21 * m22 * m23;
	}
	return m1 - m2;
  } catch (const char *str) {cerr << str << endl; return MatrixXcd::Zero(2*N,2*N);}
}
MatrixXcd Potential::triplefield_jacobian1 () const {
  int N = (Li+1) * (Li+1);
  try {
    if (N % 2 == 0)
      throw "ERROR<triplefield_jacobian> Size of Vector is not ODD.";
    if (num_f != 3)
      throw "ERROR<triplefield_jacobian> Number of superfields is not 3.";

	MatrixXcd w11 = MatrixXcd::Zero(2*N, 2*N);
	{ MatrixXcd w; w = field2matrix(0);
	  w11.block(0,0, N,N) = w; w11.block(N,N, N,N) = w.adjoint();}
	MatrixXcd w12 = MatrixXcd::Zero(2*N, 2*N);
	{ MatrixXcd w; w = field2matrix(1);
	  w12.block(0,0, N,N) = w; w12.block(N,N, N,N) = w.adjoint();}
	MatrixXcd w13 = MatrixXcd::Zero(2*N, 2*N);
	{ MatrixXcd w; w = field2matrix(2);
	  w13.block(0,0, N,N) = w; w13.block(N,N, N,N) = w.adjoint();}
	MatrixXcd w22 = MatrixXcd::Zero(2*N, 2*N);
	{ MatrixXcd w; w = field2matrix(3);
	  w22.block(0,0, N,N) = w; w22.block(N,N, N,N) = w.adjoint();}
	MatrixXcd w23     = MatrixXcd::Zero(2*N, 2*N);
	MatrixXcd w23_inv = MatrixXcd::Zero(2*N, 2*N);
	{ MatrixXcd w; w = field2matrix(4);
	  PartialPivLU<MatrixXcd> lu(w);
	  MatrixXcd w_inv; w_inv = lu.inverse();
	  w23.block(0,0, N,N)     = w;     w23.block(N,N, N,N)     = w.adjoint();
	  w23_inv.block(0,0, N,N) = w_inv; w23_inv.block(N,N, N,N) = w_inv.adjoint();}
	MatrixXcd w33 = MatrixXcd::Zero(2*N, 2*N);
	{ MatrixXcd w; w = field2matrix(5);
	  w33.block(0,0, N,N) = w; w33.block(N,N, N,N) = w.adjoint();}
	{
	  int L = Li+1;
	  MatrixXcd p = MatrixXcd::Zero(N, N);
	  for (int i = 0; i < N; ++i) {
		p(i, i) = 2.0 * M_PI * double(Li) * (IM * double(L/2 - i/L) + double(L/2 - i%L));
	  }
	  w11.block(N,0, N,N) = p; w11.block(0,N, N,N) = - p.adjoint();
	  w22.block(N,0, N,N) = p; w22.block(0,N, N,N) = - p.adjoint();
	  w33.block(N,0, N,N) = p; w33.block(0,N, N,N) = - p.adjoint();
	}

	MatrixXcd m1, m2;
	m1 = w23 - w33 * w23_inv * w22;
	{
	  MatrixXcd m2_1, m2_2, m2_3, m2_4;
	  m2_1 = w11 - w13 * w23_inv * w12;
	  m2_2 = w12 - w13 * w23_inv * w22;
	  {PartialPivLU<MatrixXcd> lu(m1); m2_3 = lu.inverse();}
	  m2_4 = w13 - w33 * w23_inv * w12;
	  m2 = m2_1 - m2_2 * m2_3 * m2_4;
	}
	return m1 * m2;
  } catch (const char *str) {cerr << str << endl; return MatrixXcd::Zero(4*N,4*N);}
}
MatrixXcd Potential::triplefield_jacobian2 () const {
  int N = (Li+1) * (Li+1);
  try {
    if (N % 2 == 0)
      throw "ERROR<triplefield_jacobian> Size of Vector is not ODD.";
    if (num_f != 3)
      throw "ERROR<triplefield_jacobian> Number of superfields is not 3.";

	MatrixXcd w11 = MatrixXcd::Zero(2*N, 2*N);
	{ MatrixXcd w; w = field2matrix(0);
	  w11.block(0,0, N,N) = w; w11.block(N,N, N,N) = w.adjoint();}
	MatrixXcd w12     = MatrixXcd::Zero(2*N, 2*N);
	MatrixXcd w12_inv = MatrixXcd::Zero(2*N, 2*N);
	{ MatrixXcd w; w = field2matrix(1);
	  PartialPivLU<MatrixXcd> lu(w);
	  MatrixXcd w_inv; w_inv = lu.inverse();
	  w12.block(0,0, N,N)     = w;     w12.block(N,N, N,N)     = w.adjoint();
	  w12_inv.block(0,0, N,N) = w_inv; w12_inv.block(N,N, N,N) = w_inv.adjoint();}
	MatrixXcd w13 = MatrixXcd::Zero(2*N, 2*N);
	{ MatrixXcd w; w = field2matrix(2);
	  w13.block(0,0, N,N) = w; w13.block(N,N, N,N) = w.adjoint();}
	MatrixXcd w22 = MatrixXcd::Zero(2*N, 2*N);
	{ MatrixXcd w; w = field2matrix(3);
	  w22.block(0,0, N,N) = w; w22.block(N,N, N,N) = w.adjoint();}
	MatrixXcd w23 = MatrixXcd::Zero(2*N, 2*N);
	{ MatrixXcd w; w = field2matrix(4);
	  w23.block(0,0, N,N) = w; w23.block(N,N, N,N) = w.adjoint();}
	MatrixXcd w33 = MatrixXcd::Zero(2*N, 2*N);
	{ MatrixXcd w; w = field2matrix(5);
	  w33.block(0,0, N,N) = w; w33.block(N,N, N,N) = w.adjoint();}
	{
	  int L = Li+1;
	  MatrixXcd p = MatrixXcd::Zero(N, N);
	  for (int i = 0; i < N; ++i) {
		p(i, i) = 2.0 * M_PI * double(Li) * (IM * double(L/2 - i/L) + double(L/2 - i%L));
	  }
	  w11.block(N,0, N,N) = p; w11.block(0,N, N,N) = - p.adjoint();
	  w22.block(N,0, N,N) = p; w22.block(0,N, N,N) = - p.adjoint();
	  w33.block(N,0, N,N) = p; w33.block(0,N, N,N) = - p.adjoint();
	}

	MatrixXcd m1 = MatrixXcd::Zero(4*N, 4*N);
	MatrixXcd m2 = MatrixXcd::Zero(4*N, 4*N);
	{
	  m1.block(0,  0, 2*N,2*N) = w12; m1.block(  0,2*N, 2*N,2*N) = w13;
	  m1.block(2*N,0, 2*N,2*N) = w23; m1.block(2*N,2*N, 2*N,2*N) = w33;
	}
	{
	  MatrixXcd m2_1 = MatrixXcd::Zero(4*N, 2*N);
	  MatrixXcd m2_2 = MatrixXcd::Zero(2*N, 4*N);
	  m2_1.block(0,0, 2*N,2*N) = w11; m2_1.block(2*N,0, 2*N,2*N) = w13;
	  m2_2.block(0,0, 2*N,2*N) = w22; m2_2.block(0,2*N, 2*N,2*N) = w23;
	  m2 = m2_1 * w12_inv * m2_2;
	}
	return m1 - m2;
  } catch (const char *str) {cerr << str << endl; return MatrixXcd::Zero(4*N,4*N);}
}
MatrixXcd Potential::multifield_jacobian () const {
  int N = (Li+1) * (Li+1);
  int Ns = 2 * N * num_f;
  try {
    if (N % 2 == 0)
      throw "ERROR<multifield_jacobian> Size of Vector is not ODD.";
    if (num_f < 1)
      throw "ERROR<multifield_jacobian> Number of superfields is not positive.";
	int L = Li+1;
	MatrixXcd jac = MatrixXcd::Zero(Ns, Ns);
	{
	  MatrixXcd iP = MatrixXcd::Zero(2*N, 2*N);
	  {
		MatrixXcd p = MatrixXcd::Zero(N, N);
		for (int i = 0; i < N/2; ++i) {
		  p(i, i) = 2.0 * M_PI * double(Li) * (IM * double(L/2 - i/L) + double(L/2 - i%L));
		  p(N-1 - i, N-1 - i) = - p(i,i);
		}
		iP.block(0,0, N,N) = p; iP.block(N,N, N,N) = - p.conjugate();
	  }
	  for (int i = 0; i < num_f; ++i) {
		int tmp = 2*N * i;
		jac.block(tmp, tmp, 2*N, 2*N) = iP;
	  }
	}

	int column_num = 0;
	for (int i = 0; i < num_f; ++i) {
	  int tmp1 = 2*N * i;
	  for (int j = i; j < num_f; ++j) {
		int tmp2 = 2*N * j;
		MatrixXcd w; w = field2matrix(column_num);
		jac.block(tmp1+N, tmp2, N,N) = w;
		jac.block(tmp1, tmp2+N, N,N) = w.adjoint();
		if (i != j) {
		  jac.block(tmp2+N, tmp1, N,N) = w;
		  jac.block(tmp2, tmp1+N, N,N) = w.adjoint();
		}
		++column_num;
	  }
	}
	return jac;
  } catch (const char *str) {cerr << str << endl; return MatrixXcd::Zero(Ns,Ns);}
}

MatrixXcd Potential::jacobian () const {
  if (num_f == 1) {return unifield_jacobian();}
  else if (num_f == 2) {return dualfield_jacobian();}
  else if (num_f == 3) {return triplefield_jacobian2();}
  return MatrixXcd::Zero( (Li+1)*(Li+1), (Li+1)*(Li+1) );
}

int Potential::sign_det (const int debug) const {
  if (num_f < 4 && num_f > 0) {
	MatrixXcd jac; jac = jacobian();
	int sign = MatrixComp::SignDet(jac, debug);
	if (num_f == 1) {return (Li%4==0)?(sign):(- sign);}
	else if (num_f == 2) {return sign;}
	else if (num_f == 3) {return sign;}
  }
  else if (num_f >=4) {
	return sign_det4multifield(debug);
  }
  return 0;
}
int Potential::sign_det4multifield (const int debug) const {
  MatrixXcd jac = multifield_jacobian();
  int ph = (num_f % 2 == 0)?1:(-1);
  return ph * MatrixComp::SignDet(jac, debug);
}

void Potential::show () const {
  cout << "# class Potential" << endl
       << "Li     = " << Li << endl
       << "num_f  = " << num_f << endl
       << "field(i, j) = " << endl
       << field << endl;
  return;
}



MatrixXd PotentialNR::nr_mat () const {
  int L = Li + 1;
  int N = L * L;

  MatrixXd p0 = MatrixXd::Zero(N, N);
  MatrixXd p1 = MatrixXd::Zero(N, N);
  for (int i = 0; i < N; ++i) {
    p0(i, i) = 2 * M_PI / double(Li) * double(L/2 - i/L);
    p1(i, i) = 2 * M_PI / double(Li) * double(L/2 - i%L);
  }

  MatrixXd NR = MatrixXd::Zero(2*N * num_f, 2*N * num_f);
  for (int i = 0; i < num_f; ++i) {
    int tmp = 2*N * i;
    //
    //  /p1 -p0\
    //  \p0  p1/
    //
    NR.block(tmp,  tmp,  N,N) =    NR.block(tmp+N,tmp+N,N,N) = p1;
    NR.block(tmp,  tmp+N,N,N) = - (NR.block(tmp+N,tmp,  N,N) = p0);
  }
  int column_num = 0;
  for (int i = 0; i < num_f; ++i) {
    int tmp1 = 2*N * i;
    for (int j = i; j < num_f; ++j) {
      int tmp2 = 2*N * j;
      MatrixXcd w2f = field2matrix(column_num) / double( Li*Li );
      MatrixXcd w2r = MatrixXcd::Zero(N, N);
      for (int l = 0; l < N; ++l) {w2r.row(l) = w2f.row(N-1 - l);}
      MatrixXd w2r_real = w2r.real();
      MatrixXd w2r_imag = w2r.imag();
      //
      //  /re  -im\
      //  \-im -re/
      //
      NR.block(tmp1,  tmp2,  N,N) += w2r_real;
      NR.block(tmp1,  tmp2+N,N,N) -= w2r_imag;
      NR.block(tmp1+N,tmp2,  N,N) -= w2r_imag;
      NR.block(tmp1+N,tmp2+N,N,N) -= w2r_real;
      if (i != j) {
		NR.block(tmp2,  tmp1,  N,N) += w2r_real;
		NR.block(tmp2,  tmp1+N,N,N) -= w2r_imag;
		NR.block(tmp2+N,tmp1,  N,N) -= w2r_imag;
		NR.block(tmp2+N,tmp1+N,N,N) -= w2r_real;
      }
      ++column_num;
    }
  }
  return NR;
}
VectorXd PotentialNR::nr_pvec () const {
  int N = (Li+1) * (Li+1);
  VectorXd NR = VectorXd::Zero(2*N * num_f);
  for (int j = 0; j < num_f; ++j) {
    int tmp = 2*N * j;
    VectorXcd w1; w1 = field2.col(j).reverse();
    VectorXcd w2; w2 = field3.col(j).reverse();
    for (int i = 0; i < N; ++i) {
      // NR(i)     =  double(k-2) * w1(i).real() + n(i).real();
      // NR(i + N) = -double(k-2) * w1(i).imag() + n(i).imag();
      NR(tmp + i)     =   w2(i).real() - w1(i).real();
      NR(tmp + i + N) = - w2(i).imag() + w1(i).imag();
    }
  }
  return NR;
}

VectorXcd PotentialNR::nrerr_pvec () const {
  int N = (Li+1) * (Li+1);
  int num_col = field2.cols();
  VectorXcd NR = VectorXcd::Zero(N * num_col);
  for (int j = 0; j < num_col; ++j) {
    int tmp = N * j;
    for (int i = 0; i < N; ++i) {NR(tmp + i) = field2(N-1 - i, j);}
  }
  return NR;
}

void PotentialNR::show () const {
  cout << "# class PotentialNR" << endl
	   << "Li     = " << Li << endl
	   << "num_f  = " << num_f << endl
	   << "field(i, j)|_{i <= j} = " << endl
	   << field << endl
	   << "field2(i) = " << endl
	   << field2 << endl
	   << "field3(i) = " << endl
	   << field3 << endl;
  return;
}

