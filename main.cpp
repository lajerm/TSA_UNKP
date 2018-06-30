
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>

#include "lapacke.h"

#include "phi4.h"
#include "subspaces.h"

#include <iomanip>

double HilbertSpace::L;
double HilbertSpace::M;
int HilbertSpace::periodicBC;

int Min(int l1, int l2)
{
    return (l1<l2)?l1:l2;
}

double _facstep(int l, int n)
{
    double val=1.0;
    for (int i=0; i <n; i++)
    {
        val *= std::sqrt(((double)l-(double)i)/((double)n-(double)i));
    }
    for (int i=n; i < l; i++)
    {
        val /= std::sqrt((double)l-(double)i);
    }
    return val;
}

double omega(int muL, int m)
{
    return std::sqrt((double)muL*(double)muL+4.0*pi*pi*(double)m*(double)m);
}

double MinCosMat(int l1, int l2, int muL, int m)
{
    double val=0.0;
    for (int n = 0; n <= Min(l1,l2); n++)
    {
        if ((l1+l2)%2==0)
        {
            val +=_facstep(l1,n)*_facstep(l2,n)*std::pow(-2.0*pi/omega(muL,m),(l1+l2)/2-n);
        } else
        {
            val +=_facstep(l1,n)*_facstep(l2,n)*std::pow(-2.0*pi/omega(muL,m),(l1+l2-1)/2-n);
        }

    }
    if ((l1+l2)%2 != 0)
    {
        val *=std::sqrt(2.0*pi/omega(muL,m));
    }
    return val;
}

int generateMxEL(int nmax, std::vector<std::vector<double> > &mx_out, int muL)
{
    for (int i =-nmax; i <= nmax; i++)  //which momentum
    {
        std::vector<double> temp((nmax+1)*(nmax+1));
        for (int j=0; j<= nmax; j++)
        {
            for (int k=0; k<=nmax; k++)
            {
                temp[j*(nmax+1)+k] = MinCosMat(j,k,muL,i);
            }
        }
        mx_out.push_back(temp);
    }
}

double CosMxEl(HilbertState bra, HilbertState ket, const std::vector<std::vector<double> > &table, int nmax)  //FIgyelem! Cos az l1+l2 páros esetben, ellenkező esetben sin!
{
	int it1 =0, it2=0;
	double num = 1.0;
	while ((it1 < bra.StateNum()) && (it2 < ket.StateNum()))
	{
		if (bra.getState(it1).k == ket.getState(it2).k)
		{
            num *= table[bra.getState(it1).k+nmax][bra.getState(it1).n*(nmax+1)+ket.getState(it2).n];
			it1++;
			it2++;
		} else
		{
			if (compare(bra.getState(it1),ket.getState(it2)))
			{
                num *= table[bra.getState(it1).k+nmax][bra.getState(it1).n*(nmax+1)+0];
				it1++;
			} else
			{
                num *= table[ket.getState(it2).k+nmax][0+ket.getState(it2).n];
				it2++;
			}
		}
	}
	if (it1==bra.StateNum())
	{
		for (int j=it2; j < ket.StateNum(); j++)
        {
            num *= table[ket.getState(j).k+nmax][0+ket.getState(j).n];
        }
	} else
	{
		for (int j=it1; j < bra.StateNum(); j++)
        {
            num *= table[bra.getState(j).k+nmax][bra.getState(j).n*(nmax+1)+0];
        }
	}
    return num;
}

int partitionVec(std::vector<int>& part_out)
{
    part_out.clear();
    part_out.push_back(1);
    for (int i = 1; i <20; i++)
    {
        int pp = 0;
        for (int j = 1; j <=i; j++)
        {
            int q = 1;
            if (j%2==1) q = 1;
            else
                q=-1;
            int gj = j*(3*j-1)/2;
            if (gj <= i)
                pp += part_out[i-gj]*q;

            gj = j*(3*j+1)/2;
            if (gj <= i)
                pp += part_out[i-gj]*q;
        }
        part_out.push_back(pp);
    }
    return 0;
}

int oddVec(const std::vector<int>& part_in, std::vector<int>& odd_out)
{
    odd_out.clear();
    odd_out.push_back(1);
     for (int i = 1; i <part_in.size(); i++)
    {
        int pp = 0;
        for (int j = 1; j <=i; j++)
        {
            int q = 1;
            if (j%2==1) q = -1;
            else
                q=1;
            int gj = j*(3*j-1);
            if (gj <= i)
                pp += part_in[i-gj]*q;

            gj = j*(3*j+1);
            if (gj <= i)
                pp += part_in[i-gj]*q;
        }
        pp += part_in[i];
        odd_out.push_back(pp);
    }
    return 0;

}

int printState(HilbertState state_in)
{
    for (int i=0; i < state_in.StateNum(); i++)
    {
        std::cout << '(' << state_in.getState(i).k << ',' << state_in.getState(i).n << ") ";
    }
    std::cout << '\n';
}

void parity_part(const HilbertSpace& HS_in, int &np, int &na, const std::vector<double> &mx_in, std::vector<std::vector<double> > &dm_out)
{
    std::vector<double> even;
    std::vector<double> odd;
    std::vector<double> evenodd;
    std::vector<int> p_index;

    int j0 = 0, j1 = 1;
    int dim = HS_in.NumStates();
    for (int i = 0; i < dim; i++)
    {

        if ((HS_in.GetState(i).ParticleNumber()%2) == 0)
        {
            p_index.push_back(j0);
            j0 += 2;
        } else
        {
            p_index.push_back(j1);
            j1 += 2;
        }
    }

    np = j0/2;
    na = j1/2;

    for (int i = 0; i < dim; i++)
    {
        for (int j=0; j < dim; j++)
        {
            if ((p_index[i]%2)==0)
            {
                if((p_index[j]%2)==0)
                {
                    even.push_back(mx_in[i*dim+j]);
                }
                else
                {
                    evenodd.push_back(mx_in[i*dim+j]);
                }
            }
            else
            {
                if((p_index[j]%2)==1)
                {
                    odd.push_back(mx_in[i*dim+j]);
                }
            }
        }
    }
    dm_out.push_back(even);
    even.clear();
    dm_out.push_back(odd);
    odd.clear();
    dm_out.push_back(evenodd);
    evenodd.clear();
    p_index.clear();
}

int matrixTrunc(std::vector<double> &mx_inout, short origWidth, short mSpaceKept)
{
    std::vector<double> temp;
    for (short i=0; i < mSpaceKept; i++)
    {
        for (short j=0; j < mSpaceKept; j++)
            temp.push_back(mx_inout[i*origWidth+j]);
    }
    mx_inout =temp;
    return 0;
}

bool compareE(const std::vector<std::vector<int> > &secinf_in)
{
    short parity = 0;
    for (int i = 0; i < 2; i++) parity += secinf_in[i][0];
    if ((parity%2)==0) return true;  // even parity sector
    return false;
}

bool compareO(const std::vector<std::vector<int> > &secinf_in)
{
    short parity = 0;
    for (int i = 0; i < 2; i++) parity += secinf_in[i][0];
    if ((parity%2)==1) return true;  // even parity sector
    return false;
}

struct genHSpace
{
    int minSpaceKept;
    int np;
    int na;

    genHSpace(int minin, int np_in, int na_in) : minSpaceKept(minin),np(np_in),na(na_in) {}
    int operator()(const std::vector<int> &subspaceVector, std::vector<std::vector<int> > &states_out)
    {
        int ind = 0;
        switch (subspaceVector[1])
        {
            case 0: ind = np; break;
            case 1: ind = na; break;
        }
        for (int i=0; i < minSpaceKept; i++)
        {
            for (int j=0; j < ind; j++)
            {
                states_out.push_back(std::vector<int> {i,j});
            }
        }
    }
};

short parity_in;
short partnum_in;

int calculateSpectrum(HilbertSpace H, double fmass, int partnum, short origWidth, short mSpaceKept, int l, std::string fname)
{
    std::fstream fs;        // file variable
    std::string s;          // will contain file name
    double dd;              // contain double read from file
    int dim=H.NumStates();  // dimension of nonminimal Hamiltonian

    // vectors to contain minimal matrices
    std::vector<double> chmp;
    std::vector<double> chmo;
    std::vector<double> shm;
    std::vector<double> evalp;
    std::vector<double> evalo;

    // Read minispace data from files
    s=fname+std::to_string(l)+"chmp";
    fs.open(s.c_str(),std::ios::in);
    while (!fs.eof()) {fs>>dd; chmp.push_back(dd);}
    fs.close();
    s=fname+std::to_string(l)+"chmo";
    fs.open(s.c_str(),std::ios::in);
    while (!fs.eof()) {fs>>dd; chmo.push_back(dd);}
    fs.close();
    s=fname+std::to_string(l)+"shm";
    fs.open(s.c_str(),std::ios::in);
    while (!fs.eof()) {fs>>dd; shm.push_back(dd);}
    fs.close();
    s=fname+std::to_string(l)+"evalp";
    fs.open(s.c_str(),std::ios::in);
    while (!fs.eof()) {fs>>dd; evalp.push_back(dd);}
    fs.close();
    s=fname+std::to_string(l)+"evalo";
    fs.open(s.c_str(),std::ios::in);
    while (!fs.eof()) {fs>>dd; evalo.push_back(dd);}
    fs.close();

    // truncate minimal matrices to dimension mSpaceKept
    matrixTrunc(chmp,origWidth,mSpaceKept);
    matrixTrunc(chmo,origWidth,mSpaceKept);
    matrixTrunc(shm,origWidth,mSpaceKept);
    // create diagonal energy matrices
    std::vector<double> mevalp(mSpaceKept*mSpaceKept);
    std::vector<double> mevalo(mSpaceKept*mSpaceKept);
    std::vector<double> midp(mSpaceKept*mSpaceKept);
    std::vector<double> mido(mSpaceKept*mSpaceKept);

    for (int i=0; i < mSpaceKept; i++)
    {
        mevalp[i*(mSpaceKept+1)] = evalp[i];
        mevalo[i*(mSpaceKept+1)] = evalo[i];
        midp[i*(mSpaceKept+1)] = 1.0;
        mido[i*(mSpaceKept+1)] = 1.0;
    }

    // generate auxiliary table of matrix elements
    std::cerr << "L=" << l << ": generateMxEl ";
    std::vector<std::vector<double> > table;
    generateMxEL(partnum, table, l);

    //generate nonminimal cos and sin matrices
    std::cerr << "OK fillStates ";
    std::vector<double> chmat(dim*dim);
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            chmat[i*dim+j]=CosMxEl(H.GetState(i),H.GetState(j),table,partnum);
        }
    }
    std::cerr << "OK AssembleMx: ";

    // generate nonminimal diagonal matrix H0
    std::vector<double> energymat(dim*dim,0.0);
    for(int i=0; i< dim; i++)
        energymat[i*(dim+1)]=H.Energy0(i);

    // separate even-odd subsectors
    std::vector<std::vector<double> > chtomb;
    std::vector<std::vector<double> > entomb;
    int np, na;
    parity_part(H, np, na, chmat, chtomb);
    parity_part(H, np, na, energymat, entomb);

    // generate nonminimal identity matrices
    std::vector<double> idp(np*np,0.0);
    std::vector<double> ido(na*na,0.0);
    for (int i=0; i< np; i++)
        idp[i*(np+1)]=1.0;
    for (int i=0; i < na; i++)
        ido[i*(na+1)]=1.0;

    // setting up data structure
    // minispace basic operators
    basic_operator_dense _mIdp(0,0,midp);
    basic_operator_dense _mIdo(1,1,mido);
    basic_operator_dense _mH0p(0,0,mevalp);
    basic_operator_dense _mH0o(1,1,mevalo);
    basic_operator_dense _mCosp(0,0,chmp);
    basic_operator_dense _mCoso(1,1,chmo);
    basic_operator_dense _mSin(0,1,shm);
    //minispace operators
    operator_dense mId {_mIdp, _mIdo};
    operator_dense mH0 {_mH0p,_mH0o};
    operator_dense mCos {_mCosp,_mCoso};
    operator_dense mSin {_mSin};

    // nonmini space basic operators
    basic_operator_dense _idp(0,0,idp);
    basic_operator_dense _ido(1,1,ido);
    basic_operator_dense _H0p(0,0,entomb[0]);
    basic_operator_dense _H0o(1,1,entomb[1]);
    basic_operator_dense _Cosp(0,0,chtomb[0]);
    basic_operator_dense _Coso(1,1,chtomb[1]);
    basic_operator_dense _Sin(0,1,chtomb[2]);
    // nonmini space operators
    operator_dense Id {_idp, _ido};
    operator_dense H0 {_H0p,_H0o};
    operator_dense Cos {_Cosp,_Coso};
    operator_dense Sin {_Sin};

    std::vector<std::vector<int> > secInf {std::vector<int> {0}, std::vector<int> {1} };
    subSpace<operator_dense> sub1(std::vector<operator_dense> {mId,mH0,mCos,mSin},secInf,std::vector<int> {mSpaceKept,mSpaceKept});
    subSpace<operator_dense> sub2(std::vector<operator_dense> {Id, H0,Cos,Sin},secInf,std::vector<int> {np,na});

    // set up full truncated Hilbert space
    std::vector<subSpace<operator_dense> > denseSubspaces{sub1,sub2};
    std::vector<subSpace<operator_sparse> > sparseSubspaces;
    StateSpace StateSp(denseSubspaces,sparseSubspaces);
    if (parity_in ==0)
    {
        StateSp.generateAllowedSubspaces(compareE);
    } else
        StateSp.generateAllowedSubspaces(compareO);

    StateSp.generateAllowedStates(genHSpace(mSpaceKept,np,na));

    // define the interaction
    std::vector<int> opmH0 {1,0};
    std::vector<int> opH0 {0,1};
    std::vector<int> opCos {2,2};
    std::vector<int> opSin {3,3};
    std::vector<int> opCosId {2,0};

    std::vector<double> coefficients {1.0,1.0,(double)l*fmass,-(double)l*fmass,-(double)l*fmass};  // to be modified
    std::vector<std::vector<int> > interaction {opmH0,opH0,opCos,opSin,opCosId};

    // compute Hamiltonian matrix
    HamiltonianDense HD(coefficients,interaction,&StateSp);
    HD.generateHamiltonian();

    std::vector<double> matrix;
    int N =0;
    HD.extractMatrix(matrix,N);

    //diagonalize
    std::cerr << "OK Diagonalize:\n ";

    double* mx = (double*)calloc(matrix.size(),sizeof(double));
    double* evals = (double*)calloc(N,sizeof(double));
    std::copy(matrix.begin(),matrix.end(),mx);
    int info;
    char cV= 'N';
    char cU= 'L';
    double worksize;
    int lwork=-1;
    LAPACK_dsyev(&cV,&cU,&N,mx,&N,evals,&worksize,&lwork,&info);
    lwork=worksize;

    double* work=(double*)calloc(lwork,sizeof(double));

    LAPACK_dsyev(&cV,&cU,&N,mx,&N,evals,work,&lwork,&info);
//    std::cout << "Last eval=" << evals[0] << '\n';
    char pc;
    if (parity_in == 0)
    {
        pc = 'e';
    } else
        pc = 'o';
    fs.open((fname+pc+std::string(".")+std::to_string(partnum_in)).c_str(),std::ios::out | std::ios::app);
    fs << l << ' ';
    for (int i=0; i < 7; i++)
        fs << std::setprecision(10) << evals[i] << ' ';
    fs << '\n';
    fs.close();
    free(work);
    free(evals);
    free(mx);

/*    fs.open("testmat",std::ios::out | std::ios::trunc);
    for (int i=0; i < mSpaceKept*(np+na); i++)
    {
        for (int j=0; j < mSpaceKept*(np+na); j++)
        {
            fs << matrix[i*mSpaceKept*(np+na)+j] << ' ';
        }
        fs << '\n';
    }
    fs.close();*/
    HD.clear();
}



int main(int argc, char* argv[])
{
    if (argc != 5)
    {
        std::cout << "usage: " << argv[0] << "<directory_in> <partnum_in> <parity_in> <mass_sign>\n";
        return 0;
    }
    std::cout << MinCosMat(0,1,1,0) << '\n';

    double M = 1.0;

    short PBC = 1;
    short partnum = atoi(argv[2]);
    partnum_in = partnum;
    double E0max = partnum*2.0*pi;
    parity_in = atoi(argv[3]);
    char* dir =argv[1];
    short origWidth = 11;
    short mSpaceKept = 7;
    short ms = atoi(argv[4]);
    std::vector<double> fmassvec {ms*0.01,ms*0.02,ms*0.04,ms*0.06,ms*0.07,ms*0.08,ms*0.09,ms*0.1,ms*0.2,ms*0.4,ms*0.7,ms*1.0};

    std::cout << "Calculating partitions: ";

    std::vector<int> partvec(1);
    std::vector<int> oddvec(1);
    partitionVec(partvec);      //the first elements of p(n)
    oddVec(partvec,oddvec);     //the first elements of a(n)

    if (PBC==-1)
        partvec = oddvec;
    for (int i = 0; i < partvec.size(); i++)
    {
        partvec[i] *= partvec[i];
        if (i > 0)
        {
            partvec[i] += partvec[i-1];
        }
        std::cout << partvec[i] << ' ';
    }

    std::cout << "Setting up basis:\n";

    HilbertSpace H;

    H.setParams(1.0,0.0,PBC);
    HilbertState vacuum;

    E0max = basis(H,vacuum,1, E0max+0.1,partvec[partnum],0,0); //partvec[partnum]

    if (E0max<0.0)
        return -1;
    std::cout << H.NumStates() << " states generated\n";
    for (int i=9; i <= 9; i++)
    {
        std::string fname = "/home/raynor/ZeroMode/"+std::string(dir)+"/m"+std::to_string(i+1)+"/th0m"+std::to_string(i+1)+"L";
        for (int l=1; l <= 20; l++)
        {
            H.setParams(l,M,PBC);
            calculateSpectrum(H, fmassvec[i],partnum,origWidth, mSpaceKept, l, fname);
        }
    }

/*    std::vector<double> mxtest;
    for (int i=0; i < H.NumStates(); i++)
    {
        for (int j=0; j < H.NumStates(); j++)
        {
            mxtest.push_back((H.GetState(i).ParticleNumber()%2)+(H.GetState(j).ParticleNumber()%2));
        }
    }
    std::vector<std::vector<double> > mxtomb;
    int np, na;
    parity_part(H, np, na, mxtest, mxtomb);
    std::cout << "np: " << np << ' ' << "na: "<< na << '\n';
    for (int i=0; i < np; i++)
    {
        for (int j=0; j < na; j++)
        {
            std::cout << mxtomb[2][i*na+j] << ' ';
        }
        std::cout << '\n';
    }*/
    return 0;
}

/*

int main()
{
    std::cout << MinCosMat(0,1,1,0) << '\n';

    double M = 1.0;
    double E0max = 11.0*2.0*pi;
    short PBC = 1;
    short partnum = 11;

    std::cout << "Calculating partitions: ";

    std::vector<int> partvec(1);
    std::vector<int> oddvec(1);
    partitionVec(partvec);      //the first elements of p(n)
    oddVec(partvec,oddvec);     //the first elements of a(n)

    if (PBC==-1)
        partvec = oddvec;
    for (int i = 0; i < partvec.size(); i++)
    {
        partvec[i] *= partvec[i];
        if (i > 0)
        {
            partvec[i] += partvec[i-1];
        }
        std::cout << partvec[i] << ' ';
    }

    std::cout << "Setting up basis:";

    HilbertSpace H;

    H.setParams(1.0,0.0,PBC);
    HilbertState vacuum;

    E0max = basis(H,vacuum,1, E0max+0.1,partvec[partnum],0,0);

    if (E0max<0.0)
        return -1;
    std::cout << H.NumStates() << " states generated\n";
    std::vector<std::vector<double> > table;
    generateMxEL(partnum, table, 1);
    printState(H.GetState(13));
    printState(H.GetState(15));
    std::cout << "Cos mx element: " << CosMxEl(H.GetState(13),H.GetState(15),table,partnum) << "\n";
    double sparsity = 0.0;
    double maxEL=0.0;
    for (int i = 0; i < partvec[partnum]; i++)
    {
        for (int j = 0; j < partvec[partnum]; j++)
        {
            double d=std::abs(CosMxEl(H.GetState(i),H.GetState(j),table,partnum));
            if (d>maxEL) maxEL =d;
        }
    }
   // sparsity /= (double)sqr(6719);
    std::cout << "Maximal element: " << maxEL;
    return 0;
}
*/

// diagonalize with LAPACK DSYEVX

/*   //diagonalize
    std::cerr << "OK Diagonalize: ";

    double* mx = (double*)calloc(matrix.size(),sizeof(double));
    double* evals = (double*)calloc(N,sizeof(double));
    int* iwork = (int*)calloc(5*N,sizeof(int));
    std::copy(matrix.begin(),matrix.end(),mx);
    int info;
    char cV= 'N';
    char cU= 'L';
    int il = 1;
    int iu = 7;
    int iM;
    int LDZ = 1;
    double abstol = 0.0;
    char range = 'I';
    double worksize;
    int lwork=-1;
    LAPACK_dsyevx(&cV,&range,&cU,&N,mx,&N,nullptr,nullptr,&il,&iu,&abstol,&iM,evals,nullptr,&LDZ,&worksize,&lwork,iwork,nullptr,&info);
    lwork=worksize;

    double* work=(double*)calloc(lwork,sizeof(double));

    LAPACK_dsyev(&cV,&cU,&N,mx,&N,evals,work,&lwork,&info);
    LAPACK_dsyevx(&cV,&range,&cU,&N,mx,&N,nullptr,nullptr,&il,&iu,&abstol,&iM,evals,nullptr,&LDZ,work,&lwork,iwork,nullptr,&info);
    std::cout << "Last eval=" << evals[0] << '\n';
    free(work);
    free(evals);
    free(iwork);
    free(mx);
*/
