/**
 * \file   field_conf.cpp
 * \brief  Execute the NR method and Compute sign of determinant
 * \see    field_class.hpp
 * \see    field_potential.hpp
 * \see    field_nicolai.hpp
 * \author Okuto Morikawa
 */
#include "field_nicolai.hpp"

/**
 * \brief Identification of scalar solutions
 * \details
 *  If the norm of difference between two Scalars < SOL_ID_MAXVAL,
 *  we regard those configurations as identical ones
 * \see   Scalar
 */
const double SOL_ID_MAXVAL = 1.0e-11;
/**
 * \brief Max error of NR iteration 
 * \see   NicolaiSol
 * \see   NicolaiSol::nr_method
 * \see   NicolaiSol::test_nr_method
 */
const double MAX_NRERR     = 1.0e-14;
/**
 * \brief   Interruption of NR iteration
 * \details Interrupt NR iteration if nr_err / min(nr_err) >= this parameter
 * \see     NicolaiSol
 * \see     NicolaiSol::nr_method
 * \see     NicolaiSol::test_nr_method
 */
const double NR_INTERRUPTION = 1.0e9;

/**
 * \brief   Structure lg_params_t
 * \param Li              Physical box size, \a N_0=N_1, which must be even integer
 * \param num_f           Number of superfields
 * \param spt             Superpotential type
 * \param testmode        test mode on/off
 * \param num_nrsolutions Number of "convergent" trials of NR method
 * \param num_nicolai     Number of configurations of Nicolai mapping
 * \param num_initialnic  Initail number of the nicolai configuration
 * \param num_loop        Maximum number of iteration
 * \param num_power       Numer of k
 * \param num_lambda      Numer of lambda
 * \param k               Power in superpotential
 * \param lambda          Coupling
 * \see     get_val
 * \see     read_input
 */
typedef struct {
  int Li;
  int num_f;
  SuperPotentialType spt;
  int testmode;
  int num_nrsolutions;
  int num_nicolai;
  int num_initialnic;
  int num_loop;
  int num_power;
  int num_coupling;
  VectorXi k;
  VectorXd lambda;
} lg_params_t;

static lg_params_t params;

/**
 * \fn    void NicolaiConf (const lg_params_t params)
 * \brief Solve a Landau--Ginzburg model by NR method;
 *        Parallel computing with openmp
 * \see lg_params_t
 */
void NicolaiConf (const lg_params_t lp) {
  try {
    if (lp.Li % 2 != 0)
      throw "ERROR<NicolaiConf> Box size is not EVEN.";

    if      (lp.spt==SuperPotentialType::AlgebraA) {
	  if (lp.num_f         != 1)
		throw "ERROR<NicolaiConf> SuperPotentialType::AlgebraA: number of fields should be 1.";
	  if (lp.k.size()      != 1)
		throw "ERROR<NicolaiConf> SuperPotentialType::AlgebraA: number of powers should be 1.";
	  if (lp.lambda.size() != 1)
		throw "ERROR<NicolaiConf> SuperPotentialType::AlgebraA: number of couplings should be 1.";
	}
    else if (lp.spt==SuperPotentialType::AlgebraD) {
	  if (lp.num_f         != 2)
		throw "ERROR<NicolaiConf> SuperPotentialType::AlgebraD: number of fields should be 2.";
	  if (lp.k.size()      != 1)
		throw "ERROR<NicolaiConf> SuperPotentialType::AlgebraD: number of powers should be 1.";
	  if (lp.lambda.size() != 2)
		throw "ERROR<NicolaiConf> SuperPotentialType::AlgebraD: number of couplings should be 2.";
	}
	else if (lp.spt==SuperPotentialType::AlgebraE) {
	  if (lp.num_f         != 2)
		throw "ERROR<NicolaiConf> SuperPotentialType::AlgebraE: number of fields should be 2.";
	  if (lp.lambda.size() != 2)
		throw "ERROR<NicolaiConf> SuperPotentialType::AlgebraE: number of couplings should be 2.";
	}
	else {
	  if (lp.num_f     < 1)
		throw "ERROR<NicolaiConf> Number of fields should be positive.";
	  if (lp.k.size() == 0)
		throw "ERROR<NicolaiConf> No values in powers.";
	  if (lp.lambda.size() == 0)
		throw "ERROR<NicolaiConf> No values in couplings";
	}

    cerr << endl << "Calc start ..." << endl;
    //--Start parallel computing (openmp)--//
    #pragma omp parallel for
    for (int in = lp.num_initialnic; in < lp.num_initialnic + lp.num_nicolai; ++in) {
      NicolaiSol nic(lp.Li, lp.num_f, in);
      double dmp = nic.nr_method(lp.k, lp.lambda, lp.num_nrsolutions, lp.num_loop, lp.spt);

      int num_sol = nic.num_sol();
      VectorXi signs(num_sol);
      for (int i = 0; i < num_sol; ++i) {signs(i) = nic.sign_det(i);}
	  if (SuperPotentialType_StdNumSol(lp.k, lp.spt) != num_sol) {
		cerr << "lm" << lp.lambda.transpose() << " k" << lp.k.transpose() << " L" << lp.Li << " n" << in
			 << "  " << "num_sol: " << num_sol << endl;
      }
      for (int i = 0; i < num_sol; ++i) {
		if (signs(i) == -1) {
		  cerr << "lm" << lp.lambda.transpose() << " k" << lp.k.transpose() << " L" << lp.Li << " n" << in
			   << "_" << i << "  " << "sign: -1" << endl;
		}
      }

      nic.nic_output(lp.k, lp.lambda, lp.num_nrsolutions, signs, dmp, nic.max_err(), lp.spt);
      nic.phi_output(in, signs);
    } //--End parallel computing (openmp)--//
    cerr << "Calc end ..." << endl << endl;
    return;
  } catch (const char *str) {cerr << str << endl; return;}
}

/**
 * \fn    void Test_SignDet (const int Li, const int num_f, const int loop)
 * \brief Test program: computation of sign_det() in Potential;
 *                      Compare signs and times,
 *                      which are computed by using two forms of jacobian
 * \param Li    Physical box size, \a N_0=N_1,
 *              which must be even integer
 * \param num_f Number of superfields
 * \param loop  Number of trials
 */
void Test_SignDet (const int Li, const int num_f, const int loop) {
  cout << "Box size        : " << Li << endl
	   << "Number of fields: " << num_f << endl
	   << "Number of loop  : " << loop << endl << endl;
  int count = 0;
  int which_fast = 0;
  VectorXd times = VectorXd::Zero(loop);
  for (int i = 0; i < loop; ++i) {
	Potential pt(Li, num_f, Distribution::Gauss_MT_unit);
	
	clock_t start1 = clock();
	int sign  = pt.sign_det();
	clock_t end1 = clock();
	double t1 = (end1 - start1);
	
	clock_t start2 = clock();
	int sign2 = pt.sign_det4multifield();
	clock_t end2 = clock();
	double t2 = (end2 - start2);
	
	if (t1 > t2) ++which_fast;
	times(i) = t1 / t2;
	if (sign != sign2) {cout << i << ": Noooooooo!!!!!!!  ("
							 << sign << ", " << sign2 << ")" << endl;}
	if (sign == -1) {
	  ++count; cout << "  #" << i << ": sign_modf  is negative" << endl;
	}
	if (sign2 == -1) {
	  cout << "  #" << i << ": sign_naive is negative" << endl;
	}
  }
  cout << endl
	   << "Positive, Negative :" << loop - count << ", " << count << endl
	   << endl
	   << "Which is fast? modf:" << loop - which_fast
	   << ", naive:" << which_fast << endl
	   << endl
	   << "Ratio of time (modf/naive): "
	   << " Average         : "
	   << times.mean() << endl
	   << " Minimum, Maximum: "
	   << times.minCoeff() << ", " << times.maxCoeff() << endl;
  return;
}


/**
 * \fn    int get_val(FILE* fp, const char *str, const char* fmt, void* val)
 */
static int get_val(FILE* fp, const char *str, const char* fmt, void* val){
  char c[128];
  if(1!=fscanf(fp,"%s",c)) {
    fprintf(stderr,"Error reading input file at %s\n",str);
    exit(1);
  }
  if(strcmp(str,c)!=0) {
    fprintf(stderr,"Error reading input file expected %s found %s\n",str,c);
    exit(1);
  }
  if(1!=fscanf(fp,fmt,val)) {
    fprintf(stderr,"Error reading input file at %s\n",str);
    fprintf(stderr,"Cannot read value format %s\n",fmt);
    exit(1);
  }
  return 0;
}

/**
 * \fn    static int read_input(char *input)
 */
static int read_input(char *input) {
  FILE* fp;
  fp=fopen(input,"r");
  if (fp==NULL) {
    fprintf(stderr,"Cannot open input file %s\n",input);
    exit(1);
  }
  int spt;
  get_val(fp, "Li",              "%i",&(params.Li));
  get_val(fp, "num_f",           "%i",&params.num_f);
  get_val(fp, "potentialtype",   "%i",&spt);
  if      (spt==0)
    params.spt = SuperPotentialType::AlgebraA;
  else if (spt==1)
    params.spt = SuperPotentialType::AlgebraD;
  else if (spt==2)
    params.spt = SuperPotentialType::AlgebraE;
  else if (spt==3)
    params.spt = SuperPotentialType::Custom;
  get_val(fp, "testmode",        "%i",&params.testmode);
  get_val(fp, "num_nrsolutions", "%i",&params.num_nrsolutions);
  get_val(fp, "num_nicolai",     "%i",&params.num_nicolai);
  get_val(fp, "num_initialnic",  "%i",&params.num_initialnic);
  get_val(fp, "num_loop",        "%i",&params.num_loop);
  get_val(fp, "num_power",       "%i",&params.num_power);
  params.k = VectorXi::Zero(params.num_power);
  for (int i = 0; i<params.num_power; ++i) {
    get_val(fp, "k",             "%i",&params.k(i));
  }
  get_val(fp, "num_coupling",    "%i",&params.num_coupling);
  params.lambda = VectorXd::Zero(params.num_coupling);
  for (int i = 0; i<params.num_coupling; ++i) {
    get_val(fp, "lambda",        "%lf",&params.lambda(i));
  }
  cerr << "Superpotential Type    : ";
  if      (params.spt==SuperPotentialType::AlgebraA)
    cerr << "Algebra A" << endl;
  else if (params.spt==SuperPotentialType::AlgebraD)
    cerr << "Algebra D" << endl;
  else if (params.spt==SuperPotentialType::AlgebraE)
    cerr << "Algebra E" << endl;
  else if (params.spt==SuperPotentialType::Custom)
    cerr << "Custom type" << endl;
  cerr << "Physical Box Size      : " << params.Li << endl
       << "Power in SuperPotential: " << params.k.transpose() << endl
       << "Coupling               : " << params.lambda.transpose() << endl << endl
       << "Configuration Number   : " << params.num_initialnic << "~" << params.num_initialnic + params.num_nicolai << endl
       << "Number of Loop         : " << params.num_loop << endl;
  return 0;
}

/**
 * \fn    int main (int argc, char* argv[])
 * \brief Main function
 */
int main (int argc, char*argv[]) {
  cout << endl << "main start ..." << endl;
  clock_t start = clock();

  // int Li = 16;
  // int num_f = 3;
  // int loop = 100;
  // Test_SignDet(Li, num_f, loop);

  if (argc != 2) {
    fprintf(stderr,"Number of arguments not correct\n");
    fprintf(stderr,"Usage: %s <infile> \n",argv[0]);
    exit(1);
  }

  read_input(argv[1]);

  if (params.testmode==1) {
    NicolaiSol nic(params.Li, params.num_f, 0);
    nic.test_nr_method(params.k,params.lambda,params.num_loop,params.num_nrsolutions,params.spt);
  }
  else
    NicolaiConf(params);

  clock_t end = clock();
  cerr << (double)(end - start) / CLOCKS_PER_SEC << "sec" << endl;
  cout << "main end ..." << endl;
}

