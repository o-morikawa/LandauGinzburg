/**
 * \file   field_class.cpp
 * \brief  Definition of the class Field
 * \see    field_class.hpp
 * \author Okuto Morikawa
 */
#include "field_class.hpp"


void Field::setgauss (const int num_col,
					  const double mean, const double dev) {
  int N = (Li+1) * (Li+1);
  field = MatrixXcd::Zero(N, num_col);
  std::random_device rnd;
  std::normal_distribution<> gauss(mean, dev);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < num_col; ++j) {
      field(i, j) = std::complex<double>(gauss(rnd), gauss(rnd));
    }
  }
  return;
}
void Field::setgaussl (const int num_col,
					   const double mean) {
  setgauss(num_col, mean, Li / M_SQRT2);
}
void Field::setgaussmt (const int num_col,
						const double mean, const double dev) {
  int N = (Li+1) * (Li+1);
  field = MatrixXcd::Zero(N, num_col);
  std::random_device rnd;
  std::mt19937_64 mt(rnd());
  std::normal_distribution<> gauss(mean, dev);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < num_col; ++j) {
      field(i, j) = std::complex<double>(gauss(mt), gauss(mt));
    }
  }
  return;
}
void Field::setgaussmtl (const int num_col,
						 const double mean) {
  setgaussmt(num_col, mean, Li / M_SQRT2);
}

void Field::setconf (const int num_col, const Distribution n,
					 const double mean, const double dev) {
  try {
    if (Li < 1 || num_f < 1 || num_col < 1)
      throw "ERROR<Field::setconf> Some parameters are not positive.";
    switch (n) {
    case Distribution::Gauss_unit:
      setgauss(num_col, mean, dev);   break;
    case Distribution::Gauss_L:
      setgaussl(num_col, mean);       break;
    case Distribution::Gauss_MT_unit:
      setgaussmt(num_col, mean, dev); break;
    case Distribution::Gauss_MT_L:
      setgaussmtl(num_col, mean);     break;
    case Distribution::Zero_conf:
      field = MatrixXcd::Zero( (Li+1)*(Li+1), num_col );
    }
  } catch (const char *str) {cerr << str << endl;
    Li = 1; num_f = 1;
    field = MatrixXcd::Zero( (Li+1)*(Li+1), 1);}
}

void Field::setconf (const VectorXcd &v, const int num_field) {
  num_f = num_field;
  int N = v.size();
  Li = (int)sqrt(N) - 1;
  try {
    if (N != (Li+1)*(Li+1))
      throw "ERROR<Field::setconf> Mismatched size of VectorXcd.";
    field = MatrixXcd::Zero(N, num_f);
    for (int i = 0; i < N; ++i) {field(i, 0) = v(i);}
  } catch (const char *str) {
    cerr << str << endl;
    ++Li; int Nnew = (Li+1)*(Li+1);
    field = MatrixXcd::Zero(Nnew, 1);
    for (int i = 0; i < N; ++i) {field(i, 0) = v(i);}
  }
}
void Field::setconf (const VectorXd &v,  const int num_field) {
  num_f = num_field;
  int Ntotal = v.size();
  int Nunit  = int(Ntotal / num_f);
  int N      = int(Nunit / 2);
  Li = (int)sqrt(N) - 1;
  MatrixXcd m(N, num_f);
  for (int j = 0; j < num_f; ++j) {
    int tmp = 2*N * j;
    for (int i = 0; i < N; ++i) {
      m(i, j) = v(tmp + i) + IM * v(tmp + i + N);
    }
  }
  try {
    if (N != (Li+1)*(Li+1))
      throw "ERROR<Field::setconf> Mismatched size of VectorXd.";
    field = m;
  } catch (const char *str) {
    cerr << str << endl;
    ++Li; int Nnew = (Li+1)*(Li+1);
    field = MatrixXcd::Zero(Nnew, num_f);
    field.block(0,0, N,num_f) = m;
  }
}
void Field::setconf (const MatrixXcd &m) {
  num_f = m.cols();
  int N = m.rows();
  Li = (int)sqrt(N) - 1;
  try {
    if (N != (Li+1)*(Li+1))
      throw "ERROR<Field::Field> Mismatched size of MatrixXcd.";
    field = m;
  } catch (const char *str) {
    cerr << str << endl;
    ++Li; int Nnew = (Li+1)*(Li+1);
    field = MatrixXcd::Zero(Nnew, num_f);
    field.block(0,0, N,num_f) = m;
  }
}
void Field::setconf (const MatrixXcd &m, const int num_field) {
  num_f = num_field;
  int num_col = m.cols();
  int N = m.rows();
  Li = (int)sqrt(N) - 1;
  try {
    if (N != (Li+1)*(Li+1))
      throw "ERROR<Field::Field> Mismatched size of MatrixXcd.";
    field = m;
  } catch (const char *str) {
    cerr << str << endl;
    ++Li; int Nnew = (Li+1)*(Li+1);
    field = MatrixXcd::Zero(Nnew, num_col);
    field.block(0,0, N,num_col) = m;
  }
}

MatrixXcd Field::vector2matrix (const int n) const {
  int L = (Li + 1);
  try {
    if (n < 0 || n >= num_f)
      throw "ERROR<Field::vector2matrix> Column number is not in [0, num_f).";
    MatrixXcd M(L, L);
    for (int i = 0; i < L; ++i) {
      for (int j = 0; j < L; ++j) {
		M(i, j) = field(L * i + j, n);
      }
    }
    return M;
  } catch (const char *str) {cerr << str << endl; return MatrixXcd::Zero(L, L);}
}
MatrixXcd Field::field2matrix (const int n) const {
  int N = (Li+1) * (Li+1);
  try {
    if (n < 0 || n >= field.cols())
      throw "ERROR<Field::field2matrix> Column number is out of range.";
    if (N % 2 == 0)
      throw "ERROR<Field::field2matrix> Size of VectorXcd is not ODD.";
    int L = (Li + 1);
    MatrixXcd M = MatrixXcd::Zero(N, N);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
		int ind_0 = L/2 + i/L - j/L;
		int ind_1 = L/2 + i%L - j%L;
		if (0 <= ind_0 && ind_0 < L && 0 <= ind_1 && ind_1 < L) {
		  M(i, j) = field(L*ind_0 + ind_1, n);
		}
      }
    }
    return M;
  } catch (const char *str) {cerr << str << endl; return MatrixXcd::Zero(N, N);}
}
MatrixXcd Field::dfield2matrix (const int n) const {
  int N = (Li+1) * (Li+1);
  try {
    if (n < 0 || n >= field.cols())
      throw "ERROR<Field::dfield2matrix> Column number is out of range.";
    if (N % 2 == 0)
      throw "ERROR<Field::dfield2matrix> Size of VectorXcd is not ODD.";
    VectorXcd v = conf(n);
    return v * v.transpose().reverse() / v(N/2);
  } catch (const char *str) {cerr << str << endl; return MatrixXcd::Zero(N, N);}
}

Field Field::conv (const int n1, const int n2) const {
  try {
    int num_col = field.cols();
    if (n1 < 0 || n1 >= num_col || n2 < 0 || n2 >= num_col)
      throw "ERROR<Field::conv> Column numbers are out of range.";
    int N = (Li+1) * (Li+1);
    if (N % 2 == 0)
      throw "ERROR<Field::conv> Size of VectorXcd is not ODD.";
    MatrixXcd mat = field2matrix(n1);
    VectorXcd cnv;
    cnv = mat * field.col(n2) / double(Li * Li);
    return Field {cnv};
  } catch (const char *str) {cerr << str << endl; return Field::Zero(Li, 1);}
}
Field Field::conjconv (const int n) const {
  try {
    int num_col = field.cols();
    if (n < 0 || n >= num_col)
      throw "ERROR<Field::conjconv> Column number is out of range.";
    int N = (Li+1) * (Li+1);
    if (N % 2 == 0)
      throw "ERROR<Field::conjconv> Size of VectorXcd is not ODD.";
    MatrixXcd m(N, 2);
    VectorXcd frev = conf(n).conjugate().reverse();
    for (int i = 0; i < N; ++i) {
      m(i, 0) = field(i, n);
      m(i, 1) = frev(i);
    }
    Field f(m);
    return f.conv(0, 1);
  } catch (const char *str) {cerr << str << endl; return Field::Zero(Li, 1);}
}
Field Field::conv_pw (const int pw, const int n) const {
  try {
    int num_col = field.cols();

    if (n < 0 || n >= num_col)
      throw "ERROR<Field::conv_pw> Column number is out of range.";
    int N = (Li+1) * (Li+1);
    if (N % 2 == 0)
      throw "ERROR<Field::conv_pw> Size of VectorXcd is not ODD.";
    if (pw < 0)
      throw "ERROR<Field::conv_pw> Number of convolutions is negative.";

    if (pw == 0) {
      int N = (Li+1) * (Li+1);
      VectorXcd v = VectorXcd::Zero(N); v(N/2) = Li * Li;
      return Field {v};
    } else if (pw == 1) {
      VectorXcd v = conf(n); return Field {v};
    } else if (pw == 2) { return conv(n, n); }

    int L = (Li + 1);
    MatrixXcd M, Mp;
    M  = field2matrix(n);
    Mp = MatrixXcd::Identity(N, N);
    for (int i = 0; i < pw-2; ++i) { Mp = M * Mp; }
    int Ln = 1;
    for (int i = 0; i < pw-1; ++i) { Ln *= Li * Li; }
    VectorXcd cnv;
    cnv = M * Mp * field.col(n) / Ln;
    return Field {cnv};
  } catch (const char *str) {cerr << str << endl; return Field::Zero(Li, 1);}
}

Field Field::combine_with(const Field &f) {
  try {
	if (field.rows()==0) {*this = f; return *this;}
	if (*this != f)
	  throw "ERROR<Field::combine_with> Sizes of Fields are not same .";
	int num_col1 = field.cols();
	int num_col2 = f.field.cols();
	int num_col_total = num_col1 + num_col2;
	int N = field.rows();
	MatrixXcd m(N, num_col_total);
	m.block(0,0,        N, num_col1) = field;
	m.block(0,num_col1, N, num_col2) = f.field;
    *this = Field {m, num_f}; return *this;
  }
  catch (const char *str) {cerr << str << endl; return *this;}
}

void Field::show () const {
  cout << "# class Field" << endl
       << "Li     = " << Li << endl
       << "field = " << endl
       << field << endl;
  return;
}

