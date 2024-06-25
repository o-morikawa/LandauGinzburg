/**
 * \file   field_nicolai.cpp
 * \brief  Definition of the classes associated with Nicolai map
 * \see    field_nicolai.hpp
 * \author Okuto Morikawa
 */
#include "field_nicolai.hpp"

VectorXd Nicolai::nr_nvec () const {
  int N = (Li+1) * (Li+1);
  VectorXd NR = VectorXd::Zero(2*N * num_f);
  for (int j = 0; j < num_f; ++j) {
    int tmp = 2*N * j;
    for (int i = 0; i < N; ++i) {
      NR(tmp + i)     = field(i, j).real();
      NR(tmp + i + N) = field(i, j).imag();
    }
  }
  return NR;
}

VectorXcd Nicolai::nrerr_nvec () const {
  int N = (Li+1) * (Li+1);
  int num_col = field.cols();
  VectorXcd NR = VectorXcd::Zero(N * num_col);
  for (int j = 0; j < num_col; ++j) {
    int tmp = N * j;
    for (int i = 0; i < N; ++i) {NR(tmp + i) = field(i, j);}
  }
  return NR;
}

void Nicolai::nic_output (const VectorXi k, const VectorXd lambda,
						  const int num_nrsol,
						  const VectorXi signs,
						  const int dmp,
						  const double max_err,
						  const SuperPotentialType spt) const {
  std::string fnm;
  if (lambda.size() == 1) {
	fnm = "lm"; fnm += std::to_string(lambda(0));
  } else {
	for(int i = 0; i < lambda.size(); ++i) {
	  fnm += "lm"; fnm += std::to_string(i);
	  fnm += "-";  fnm += std::to_string(lambda(i));
	}
  }
  if (k.size() == 1) {
	fnm += "k"; fnm += std::to_string(k(0));
  } else {
	for (int i = 0; i < k.size(); ++i) {
	  fnm += "k"; fnm += std::to_string(i);
	  fnm += "-"; fnm += std::to_string(k(i));
	}
  }
  fnm += "L"; fnm += std::to_string(Li);
  fnm += "n"; fnm += std::to_string(in);

  std::ofstream fo_n; std::string fnm_n;
  fnm_n = "confs/nicolai_";
  if      (spt==SuperPotentialType::AlgebraA)
    fnm_n += "spt-a_";
  else if (spt==SuperPotentialType::AlgebraD)
    fnm_n += "spt-d_";
  else if (spt==SuperPotentialType::AlgebraE)
    fnm_n += "spt-e_";
  else if (spt==SuperPotentialType::Custom)
    fnm_n += "spt-c_";
  fnm_n += fnm; fnm_n += ".dat";
  fo_n.open(fnm_n);
  
  fo_n << "# Superpotential Type" << endl;
  if      (spt==SuperPotentialType::AlgebraA)
    fo_n << "AlgebraA" << endl;
  else if (spt==SuperPotentialType::AlgebraD)
    fo_n << "AlgebraD" << endl;
  else if (spt==SuperPotentialType::AlgebraE)
    fo_n << "AlgebraE" << endl;
  else if (spt==SuperPotentialType::Custom)
    fo_n << "Custom" << endl;
  
  fo_n << "# Physical Box Size, Power, Coupling" << endl
       << Li << ", " << k.transpose() << ", " << lambda.transpose() << endl
       << "# Nicolai Configuration Number" << endl
       << in << endl
       << "# Number of Solutions and Dumped-Solutions" << endl
       << signs.size() << ", "
       << dmp << " / " << (num_nrsol + dmp) << endl
       << "# Signs of Determinants" << endl
       << signs.transpose() << endl
       << "# Max norm of residue" << endl
       << max_err << endl
       << "# Configuration of Gaussian random numbers" << endl
       << std::setprecision(20) << field;
  fo_n.close();
  return;
}

void Nicolai::show () const {
  cout << "# class Nicolai" << endl
       << "Li     = " << Li << endl
       << "in     = " << in << endl
       << "field = " << endl
       << field << endl;
  return;
}


bool Scalar::is_identical (const Scalar &f) const {
  double err = (field - f.field).norm() / field.norm();
  if (err < sol_id_maxval) {return true;}
  else {return false;}
}

VectorXcd Scalar::nrerr_svec () const {
  int L = Li + 1;
  int N = L * L;
  int num_col = field.cols();
  VectorXcd p = VectorXcd::Zero(N);
  for (int i = 0; i < N; ++i) {
    p(i) = 2 * M_PI / double(Li)
      * (IM * double(L/2 - i/L) + double(L/2 - i%L));
  }
  VectorXcd NR = VectorXcd::Zero(N * num_col);
  for (int j = 0; j < num_col; ++j) {
    int tmp = N * j;
    for (int i = 0; i < N; ++i) {NR(tmp + i) = p(i) * field(i, j);}
  }
  return NR;
}


Potential Scalar::superpotential_typealgebraA () const {
  try {
	if (num_f != 1)
	  throw "ERROR<Scalar::superpotential> Case SuperPotentialType::AlgebraA: number of field should be 1.";
	if (k.size() != 1)
	  throw "ERROR<Scalar::superpotential> Case SuperPotentialType::AlgebraA: size of VectorXi for powers should be 1.";
	if (lambda.size() != 1)
	  throw "ERROR<Scalar::superpotential> Case SuperPotentialType::AlgebraA: size of VectorXd for couplings should be 1.";
	//
	Field f = conv_pw(k(0)-2); f *= double(k(0)-1) * lambda(0);
	return Potential {f.conf(), num_f};
	//
  } catch (const char *str) {
	cerr << str << endl; return Potential {Li, num_f};}
}
Potential Scalar::superpotential_typealgebraD () const {
	try {
	  if (num_f != 2)
		throw "ERROR<Scalar::superpotential> Case SuperPotentialType::AlgebraD: number of fields should be 2.";
	  if (k.size() != 1)
		throw "ERROR<Scalar::superpotential> Case SuperPotentialType::AlgebraD: size of VectorXi for powers should be 1.";
	  if (lambda.size() != 2)
		throw "ERROR<Scalar::superpotential> Case SuperPotentialType::AlgebraD: size of VectorXd for couplings should be 2.";
	  //
	  Field Xk2 = conv_pw(k(0)-2,0); Xk2 *= double(k(0)-1) * lambda(0);
	  Field X1  = conv_pw(1,     0); X1  *= lambda(1);
	  Field Y1  = conv_pw(1,     1); Y1  *= lambda(1);
	  Field res = Xk2;
	  res.combine_with(Y1); res.combine_with(X1);
	  return Potential {res.conf(), num_f};
	  //
	} catch (const char *str) {
	  cerr << str << endl; return Potential {Li, num_f};}
}
Potential Scalar::superpotential_typealgebraE () const {
	try {
	  if (num_f != 2)
		throw "ERROR<Scalar::superpotential> Case SuperPotentialType::AlgebraE (E_7): number of fields should be 2.";
	  if (lambda.size() != 2)
		throw "ERROR<Scalar::superpotential> Case SuperPotentialType::AlgebraE (E_7): size of VectorXd for couplings should be 2.";
	  //
	  Field X1 = conv_pw(1, 0); Field Y1 = conv_pw(1, 1);
	  //Field XY;
	  //{Field tmp = X1; tmp.combine_with(Y1); XY = tmp.conv(0,1);}
	  Field XY = conv(0,1);
	  Field Y2 = conv_pw(2, 1);
	  X1 *= 2.0 * lambda(0); Y2 *= lambda(1); XY *= 2.0 * lambda(1);
	  Field res = X1;
	  res.combine_with(Y2); res.combine_with(XY);
	  return Potential {res.conf(), num_f};
	  //
	} catch (const char *str) {
	  cerr << str << endl; return Potential {Li, num_f};}
}

Potential Scalar::superpotential () const {
  switch (spt_type){
  default:
  case SuperPotentialType::AlgebraA:
	return superpotential_typealgebraA();
	break;
  case SuperPotentialType::AlgebraD:
	return superpotential_typealgebraD();
	break;
  case SuperPotentialType::AlgebraE:
	return superpotential_typealgebraE();
	break;
  case SuperPotentialType::Custom:
	return superpotential_typecustom();
	break;
  }
}

PotentialNR Scalar::superpotential_nr_typealgebraA () const {
	try {
	  if (num_f != 1)
		throw "ERROR<Scalar::superpotential_nr> Case SuperPotentialType::AlgebraA: number of field should be 1.";
	  if (k.size() != 1)
		throw "ERROR<Scalar::superpotentia_nrl> Case SuperPotentialType::AlgebraA: size of VectorXi for powers should be 1.";
	  if (lambda.size() != 1)
		throw "ERROR<Scalar::superpotential_nr> Case SuperPotentialType::AlgebraA: size of VectorXd for couplings should be 1.";
	  //
	  Field f1 = conv_pw(k(0)-2); f1 *= lambda(0);
	  Field f2;
	  {Field ff2 = f1; ff2.combine_with(*this); f2 = ff2.conv(0,1);}
	  f1 *= double(k(0)-1);
	  MatrixXcd m2 = f2.conf();
	  return PotentialNR {num_f, f1.conf(), m2, double(k(0)-1) * m2};
	  //
	} catch (const char *str) {
	  cerr << str << endl; return PotentialNR {Li, num_f};}
}
PotentialNR Scalar::superpotential_nr_typealgebraD () const {
	try {
	  if (num_f != 2)
		throw "ERROR<Scalar::superpotential_nr> Case SuperPotentialType::AlgebraD: number of fields should be 2.";
	  if (k.size() != 1)
		throw "ERROR<Scalar::superpotential_nr> Case SuperPotentialType::AlgebraD: size of VectorXi for powers should be 1.";
	  if (lambda.size() != 2)
		throw "ERROR<Scalar::superpotential_nr> Case SuperPotentialType::AlgebraD: size of VectorXd for couplings should be 2.";
	  //
	  Field X1  = conv_pw(1,     0);
	  Field Y1  = conv_pw(1,     1);
	  Field Y2  = conv_pw(2,     1);
	  Field Xk2 = conv_pw(k(0)-2,0);
	  Field Xk1;
	  {Field tmp = X1; tmp.combine_with(Xk2); Xk1 = tmp.conv(0,1);}
	  //Field XY;
	  //{Field tmp = X1; tmp.combine_with(Y1); XY = tmp.conv(0,1);}
	  Field XY = conv(0,1);
	  // res1: (k0 - 1) lambda0 X^(k0-2)
	  //       lambda1 Y
	  //       lambda1 X
	  X1  *= lambda(1); Y1  *= lambda(1);
	  Xk2 *= double(k(0)-1) * lambda(0);
	  Field res1 = Xk2;
	  res1.combine_with(Y1); res1.combine_with(X1);
	  // res2: lambda0 X^(k0-1) + lambda1 Y^2 / 2
	  //       lambda1 XY
	  Y2  *= lambda(1) / 2.0; XY  *= lambda(1); Xk1 *= lambda(0);
	  Field res2;
	  {MatrixXcd m2 = Xk1.conf() + Y2.conf(); res2 = Field {m2};}
	  res2.combine_with(XY);
	  // res3: (k0-1) lambda0 X^(k0-1) + lambda1 Y^2
	  //       2 lambda1 XY
	  Y2  *= 2.0; XY  *= 2.0; Xk1 *= double(k(0)-1);
	  Field res3;
	  {MatrixXcd m3 = Xk1.conf() + Y2.conf(); res3 = Field {m3};}
	  res3.combine_with(XY);
	  return PotentialNR {num_f, res1.conf(), res2.conf(), res3.conf()};
	  //
	} catch (const char *str) {
	  cerr << str << endl; return PotentialNR {Li, num_f};}
}
PotentialNR Scalar::superpotential_nr_typealgebraE () const {
	try {
	  if (num_f != 2)
		throw "ERROR<Scalar::superpotential_nr> Case SuperPotentialType::AlgebraE (E_7): number of fields should be 2.";
	  if (lambda.size() != 2)
		throw "ERROR<Scalar::superpotential_nr> Case SuperPotentialType::AlgebraE (E_7): size of VectorXd for couplings should be 2.";
	  //
	  Field X1 = conv_pw(1,0);
	  Field X2 = conv_pw(2,0);
	  Field Y1 = conv_pw(1,1);
	  Field Y2 = conv_pw(2,1);
	  Field Y3;
	  {Field tmp = Y1; tmp.combine_with(Y2); Y3 = tmp.conv(0,1);}
	  //Field XY;
	  //{Field tmp = X1; tmp.combine_with(Y1); XY = tmp.conv(0,1);}
	  Field XY = conv(0,1);
	  Field XYY;
	  {Field tmp = X1; tmp.combine_with(Y2); XYY = tmp.conv(0,1);}
	  Field YXY;
	  {Field tmp = Y1; tmp.combine_with(XY); YXY = tmp.conv(0,1);}
	  Field XY2;
	  {
		MatrixXcd tmp = ( XYY.conf() + 2.0 * YXY.conf() ) / 3.0;
		XY2 = Field {tmp};
	  }
	  // res1: 2 lambda0 X
	  //       lambda1 Y^2
	  //       2 lambda1 XY
	  X1 *= 2.0 * lambda(0); Y2 *= lambda(1); XY *= 2.0 * lambda(1);
	  Field res1 = X1;
	  res1.combine_with(Y2); res1.combine_with(XY);
	  // res2: lambda0 X^2 + lambda1 Y^3 / 3
	  //       lambda1 X Y^2
	  X2 *= lambda(0); Y3 *= lambda(1) / 3.0; XY2 *= lambda(1);
	  Field res2;
	  {MatrixXcd m2 = X2.conf() + Y3.conf(); res2 = Field {m2};}
	  res2.combine_with(XY2);
	  // res3: 2 lambda0 X^2 + lambda1 Y^3
	  //       3 lambda1 X Y^2
	  X2 *= 2.0; Y3 *= 3.0; XY2 *= 3.0;
	  Field res3;
	  {MatrixXcd m3 = X2.conf() + Y3.conf(); res3 = Field {m3};}
	  res3.combine_with(XY2);
	  return PotentialNR {num_f, res1.conf(), res2.conf(), res3.conf()};
	  //
	} catch (const char *str) {
	  cerr << str << endl; return PotentialNR {Li, num_f};}
}

PotentialNR Scalar::superpotential_nr () const {
  switch (spt_type){
  default:
  case SuperPotentialType::AlgebraA:
	return superpotential_nr_typealgebraA();
	break;
  case SuperPotentialType::AlgebraD:
	return superpotential_nr_typealgebraD();
	break;
  case SuperPotentialType::AlgebraE:
	return superpotential_nr_typealgebraE();
	break;
  case SuperPotentialType::Custom:
	return superpotential_nr_typecustom();
	break;
  }
}


Scalar Scalar::nr_loop (const Nicolai &nic) {
  PotentialNR p_nr = superpotential_nr();
  VectorXd v;
  v = p_nr.nr_pvec() + nic.nr_nvec();
  MatrixXd m = p_nr.nr_mat();
  VectorXd v_res = m.partialPivLu().solve(v);
  *this = Scalar {*this, v_res}; return *this;
}

VectorXcd Scalar::nr_error_vec (const Nicolai &nic) const {
  PotentialNR ptn = superpotential_nr();
  VectorXcd v = nrerr_svec(); // v = (p * v) //
  VectorXcd w = ptn.nrerr_pvec();
  VectorXcd n = nic.nrerr_nvec();
  return (v + w.conjugate() - n) / n.norm();
};

void Scalar::scl_output (const int in, const int num_sol, const int i,
						 const int sign, const double err) const {
  std::string fnm;
  if (lambda.size() == 1) {
	fnm = "lm"; fnm += std::to_string(lambda(0));
  } else {
	for(int i = 0; i < lambda.size(); ++i) {
	  fnm += "lm"; fnm += std::to_string(i);
	  fnm += "-";  fnm += std::to_string(lambda(i));
	}
  }
  if (k.size() == 1) {
	fnm += "k"; fnm += std::to_string(k(0));
  } else {
	for (int i = 0; i < k.size(); ++i) {
	  fnm += "k"; fnm += std::to_string(i);
	  fnm += "-"; fnm += std::to_string(k(i));
	}
  }
  fnm += "L"; fnm += std::to_string(Li);
  fnm += "n"; fnm += std::to_string(in);

  std::ofstream fo_p; std::string fnm_p;
  fnm_p = "confs/phi_";
  if      (spt_type==SuperPotentialType::AlgebraA)
    fnm_p += "spt-a_";
  else if (spt_type==SuperPotentialType::AlgebraD)
    fnm_p += "spt-d_";
  else if (spt_type==SuperPotentialType::AlgebraE)
    fnm_p += "spt-e_";
  else if (spt_type==SuperPotentialType::Custom)
    fnm_p += "spt-c_";
  fnm_p += fnm;
  fnm_p += "_"; fnm_p += std::to_string(i); fnm_p += ".dat";
  fo_p.open(fnm_p);
  
  fo_p << "# Superpotential Type" << endl;
  if      (spt_type==SuperPotentialType::AlgebraA)
    fo_p << "AlgebraA" << endl;
  else if (spt_type==SuperPotentialType::AlgebraD)
    fo_p << "AlgebraD" << endl;
  else if (spt_type==SuperPotentialType::AlgebraE)
    fo_p << "AlgebraE" << endl;
  else if (spt_type==SuperPotentialType::Custom)
    fo_p << "Custom" << endl;
  
  fo_p << "# Physical Box Size, Power, Coupling" << endl
       << Li << ", " << k.transpose() << ", " << lambda.transpose() << endl
       << "# Nicolai Configuration Number" << endl
       << in << endl
       << "# Solution Configuration Number" << endl
       << i << " / " << num_sol << endl
       << "# Sign of Determinant" << endl << sign << endl
       << "# Norm of residue" << endl << err << endl 
       << "# Configuration of a scalar field solution" << endl
       << std::setprecision(20) << field;
  fo_p.close();
}

void Scalar::show () const{
  cout << "# class Scalar" << endl
       << "Li     = " << Li << endl
       << "k      = " << k.transpose()  << endl
       << "lambda = " << lambda.transpose() << endl
       << "field = " << endl
       << field << endl;
  return;
}


NicolaiSol &NicolaiSol::add_sol(const Scalar &f, const double &err) {
  try {
    if (fs.size() == 0) {
      fs.push_back(f); errs.push_back(err); return *this;
    }
    if (fs[0] != f)
      throw "ERROR<NicolaiSol::add_sol> Sizes of Scalars are not same.";
    bool id = true;
    for(int i = 0; i < fs.size(); ++i) {
      if ( fs[i].is_identical(f) ) {id = false; break;}
    }
    if (id) {fs.push_back(f); errs.push_back(err);} return *this;
  } catch (const char *str) {cerr << str << endl; return *this;}
};
NicolaiSol &NicolaiSol::add_sol(const NicolaiSol &sol) {
  for (int i = 0; i < sol.fs.size(); ++i) {
    add_sol(sol.fs[i], sol.errs[i]);
  }
  return *this;
};

double NicolaiSol::nr_method (const VectorXi k, const VectorXd lambda,
							  const int num_nrsol, const int LOOP,
							  const SuperPotentialType t) {
  try {
	if (Li % 2 != 0)
	  throw "ERROR<NicolaiSol::nr_method> Box size is not EVEN.";
	int dmp = 0;
	for (int i = 0; i < num_nrsol; ++i) {
	  Scalar scl(Li, num_f, k, lambda, t);
	  double err; double min_err;
	  for (int j = 0; j < LOOP; ++j) {
		scl.nr_loop(*this); err = scl.nr_error(*this);
		if (err < MAX_NRERR) break;
		if (j == 0 || err < min_err) min_err = err;
		else if (err >= min_err * nr_interruption) break;
	  }
	  if (err < MAX_NRERR) add_sol(scl, err);
	  else {++dmp; --i;}
	}
	return dmp;
  } catch (const char *str) {cerr << str << endl; return -1;}
}

double NicolaiSol::nr_method_core (Scalar &scli,
							  const VectorXi k, const VectorXd lambda,
							  const int LOOP,
							  const SuperPotentialType t) {
  try {
	if (Li % 2 != 0)
	  throw "ERROR<NicolaiSol::nr_method_core> Box size is not EVEN.";
	Scalar scl;
	int dmp = -1;
	double err; double min_err;
	do {
		++dmp;
		scl = Scalar(Li, num_f, k, lambda, t);
		for (int j = 0; j < LOOP; ++j) {
			scl.nr_loop(*this); err = scl.nr_error(*this);
			if (err < MAX_NRERR) break;
			if (j == 0 || err < min_err) min_err = err;
			else if (err >= min_err * nr_interruption) break;
		}
	} while (err >= MAX_NRERR);
	scli = scl;
	return dmp;
  } catch (const char *str) {cerr << str << endl; return -1;}
}

double NicolaiSol::nr_method_omp (const VectorXi k, const VectorXd lambda,
							  const int num_nrsol, const int LOOP,
							  const SuperPotentialType t) {
  try {
	if (Li % 2 != 0)
	  throw "ERROR<NicolaiSol::nr_method_omp> Box size is not EVEN.";
	int dmp = 0;
	vector<Scalar> scls(num_nrsol);
	vector<int> dmps(num_nrsol);
	//--Start parallel computing (openmp)--//
	#pragma omp parallel for
	for (int i = 0; i < num_nrsol; ++i) {
	  dmps[i] = nr_method_core(scls[i], k, lambda, LOOP, t);
	} //--End parallel computing (openmp)--//
	for (int i = 0; i < num_nrsol; ++i) {
	  add_sol(scls[i], scls[i].nr_error(*this));
	  dmp = dmp + dmps[i];
	}
	return dmp;
  } catch (const char *str) {cerr << str << endl; return -1;}
}

void NicolaiSol::test_nr_method (const VectorXi k, const VectorXd lambda,
								 const int LOOP, const int TRIAL,
								 const SuperPotentialType t) {
  try {
	if (Li % 2 != 0)
	  throw "ERROR<NicolaiSol::test_nr_method> Box size is not EVEN.";
	double err;
	for (int j = 0; j < TRIAL; ++j) {
	  cout << endl
		   << "##########################################" << endl
		   << j+1 << "-th trial" << endl
		   << "##########################################" << endl
		   << "Physical Box Size      : " << Li << endl
		   << "Power in SuperPotential: " << k.transpose()  << endl
		   << "Coupling               : " << lambda.transpose()
		   << endl << endl;

	  Scalar scl(Li, num_f, k, lambda, t);
	  err = scl.nr_error(*this);
	  cout << "Initial Error of NR " << endl
		   << err << endl << endl;

	  clock_t start_nr = clock();

	  double min_err;
	  for (int i = 0; i < LOOP; ++i) {
		cout << "NR iteration " << std::setw(3) << i+1 << ":  ";
		Scalar tmp = scl;
		scl.nr_loop(*this); err = scl.nr_error(*this);
		int tmp_id = (scl.is_identical(tmp))?(1):(0);
		cout << "Error of NR "
			 << std::setprecision(15) << err
			 << " (" << tmp_id << ")" << endl;
		if (err < MAX_NRERR) break;
		if (i == 0 || err < min_err) min_err = err;
		else if (err >= min_err * nr_interruption) break;
	  }
	  if (err < MAX_NRERR) add_sol(scl, err);

	  clock_t end_nr = clock();
	  cout << (double)(end_nr-start_nr) / CLOCKS_PER_SEC << "sec"
		   << endl << endl;

	  if (err < MAX_NRERR) break;
	}
	return;
  } catch (const char *str) {cerr << str << endl; return;}
}

void NicolaiSol::phi_output (const int in, const VectorXi signs) const {
  int num_sol = fs.size();
  for (int i = 0; i < num_sol; ++i) {
    fs[i].scl_output(in, num_sol, i, signs(i), errs[i]);
  }
}

void NicolaiSol::show() const {
  try {
    int fssize = fs.size();
    if (fssize == 0)
      throw "ERROR<NicolaiSol::show> No entries in this object.";
    cout << "# class NicolaiSol with " << fssize << " Scalars" << endl;
    for (int i = 0; i < fssize; ++i) {
      cout << "# entry: " << i << endl << "#"; fs[i].show();
    }
    return;
  } catch (const char *str) {cerr << str << endl; return;}
}

int SuperPotentialType_StdNumSol (const VectorXi k, const SuperPotentialType t) {
  if (t==SuperPotentialType::AlgebraA)
	return k(0)-1;
  else if (t==SuperPotentialType::AlgebraD)
	return k(0)+1;
  else if (t==SuperPotentialType::AlgebraE)
	return 7;
  else if (t==SuperPotentialType::Custom)
	return SuperPotentialTypeCustom_StdNumSol;
}
