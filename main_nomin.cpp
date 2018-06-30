
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>

#include "lapacke.h"

#include "phi4.h"

#include <iomanip>

#include "coo.h"
#include "utils.h"
#include "PRIMME_int.h"
#include "phi_n.h"


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


short parity_in;
short partnum_in;

int dense2sparse(int k, const std::vector<double> &d_in, sparse_matrix &sp_out)  //k: number of columns (width)
{
    COO coo;
    for (int i = 0; i< d_in.size(); i++ )
    {
        if (d_in[i]!= 0.0)
        {
            coo.row_ind = i/k;
            coo.col_ind = i%k;
            coo.val=d_in[i];
            sp_out.PushBack(coo);
        }
    }
    return 0;
}

int calculateSpectrumTh(HilbertSpace He, HilbertSpace Ho, double fmass, double th, int partnum, int l, std::string fname)
{
    std::fstream fs;        // file variable
    std::string s;          // will contain file name
    double dd;              // contain double read from file
    int dime=He.NumStates();
    int dimo=Ho.NumStates();
    int dim = dime+dimo;
   // generate auxiliary table of matrix elements
    std::cerr << "L=" << l << ": generateMxEl ";

    std::vector<double> zmu(20);
    fs.open("ztest",std::ios::in);
    for (int i=1; i<= 20; i++)
        fs >> zmu[i-1];
    fs.close();

    int maxnum = 0;
    //calculate maximal occupation number
    for (int i=0; i < He.NumStates(); i++)
    {
        for (int j =0; j < He.GetState(i).StateNum(); j++)
        {
            if ((He.GetState(i).getState(j).k==0)&&(He.GetState(i).getState(j).n>maxnum))
                maxnum = He.GetState(i).getState(j).n;
        }
    }

    for (int i=0; i < Ho.NumStates(); i++)
    {
        for (int j =0; j < Ho.GetState(i).StateNum(); j++)
        {
            if ((Ho.GetState(i).getState(j).k==0)&&(Ho.GetState(i).getState(j).n>maxnum))
                maxnum = Ho.GetState(i).getState(j).n;
        }
    }

    if (partnum > maxnum) maxnum = partnum;
    std::cout << "maxnum: " << maxnum << '\n';
    maxnum +=1;

    // generate auxiliary matrix element table
    std::vector<std::vector<double> > table;
    generateMxEL(maxnum, table, l);

    std::cerr << "OK fillStates ";
    sparse_matrix chmat;
    COO coo;
    for (int i = 0; i < dime; i++)
    {
        coo.row_ind = i;
        for (int j = 0; j < dime; j++)
        {
            coo.col_ind = j;

            coo.val=zmu[l-1]*fmass*l*CosMxEl(He.GetState(i),He.GetState(j),table,maxnum);
            if (i==j) coo.val += He.Energy0(i);
            chmat.PushBack(coo);
        }
    }
    for (int i=0; i < dimo; i++)
    {
        coo.row_ind = dime+i;
        for (int j = 0; j < dimo; j++)
        {
            coo.col_ind = dime+j;
            coo.val=zmu[l-1]*fmass*l*CosMxEl(Ho.GetState(i),Ho.GetState(j),table,maxnum);
            if (i==j) coo.val += Ho.Energy0(i);
            chmat.PushBack(coo);
        }
    }
    COO coo2;

    interact_params ip;
    ip.coefficient=std::sqrt(1/(4.0*pi))*th;
    ip.exponent=1;
    phi_n_interaction phiint = phi_n_interaction(ip,&He);

    for (int i=0; i < dime; i++)
    {
        coo.row_ind = i;
        coo2.col_ind = i;
        for (int j=0; j < dimo; j++)
        {
            coo.col_ind = dime+j;
            coo2.row_ind = dime+j;
            coo.val = phiint.MatrixElement(He.GetState(i),Ho.GetState(j));
            coo2.val = coo.val;
            chmat.PushBack(coo);
            chmat.PushBack(coo2);
        }

    }
    chmat.sort();
    std::cerr << "OK Diagonalize:\n ";

    double* evals = (double*)calloc(dim,sizeof(double));
    //PRIMME diagonalization

    double* evecs =(double*)calloc(dim*10,sizeof(double));
    comp_evals(dim,0.1, 7, chmat, evals, evecs);
    free(evecs);

    fs.open((fname+std::string(".")+std::to_string(partnum_in)).c_str(),std::ios::out | std::ios::app);
    fs << l << ' ';
    for (int i=0; i < 7; i++)
        fs << std::setprecision(10) << evals[i] << ' ';
    fs << '\n';
    fs.close();

    free(evals);
}

int calculateSpectrum(HilbertSpace H, double fmass, int partnum, int l, std::string fname)
{
    std::fstream fs;        // file variable
    std::string s;          // will contain file name
    double dd;              // contain double read from file
    int dim=H.NumStates();  // dimension of nonminimal Hamiltonian


    // generate auxiliary table of matrix elements
    std::cerr << "L=" << l << ": generateMxEl ";

    std::vector<double> zmu(20);
    fs.open("ztest",std::ios::in);
    for (int i=1; i<= 20; i++)
        fs >> zmu[i-1];
    fs.close();
    for (int i=21; i <=40; i++) zmu.push_back(1.0);
    // calculate highest zeromode occupation number
    int maxnum = 0;
    for (int i=0; i < H.NumStates(); i++)
    {
        for (int j =0; j < H.GetState(i).StateNum(); j++)
        {
            if ((H.GetState(i).getState(j).k==0)&&(H.GetState(i).getState(j).n>maxnum))
                maxnum = H.GetState(i).getState(j).n;
        }
    }
    if (partnum > maxnum) maxnum = partnum;
    std::cout << "maxnum: " << maxnum << '\n';
    maxnum +=1;
    std::vector<std::vector<double> > table;
    generateMxEL(maxnum, table, l);
/*    table[partnum].resize(257*257);
    std::fstream fs1;
    std::fstream fs2;
    std::fstream fs3;
    fs1.open((std::string("zm/l")+std::to_string(l)+std::string("mee")).c_str(),std::ios::in);
    fs2.open((std::string("zm/l")+std::to_string(l)+std::string("moo")).c_str(),std::ios::in);
    fs3.open((std::string("zm/l")+std::to_string(l)+std::string("meo")).c_str(),std::ios::in);
    for (int i=0; i <= 128; i++)
        for (int j=0; j <= 128; j++)
        {
            fs1 >> table[partnum][(2*i)*257+(2*j)];
            if (j>0)
            {
                fs3 >> table[partnum][(2*i)*257+(2*j-1)];
                table[partnum][(2*j-1)*257+(2*i)] = table[partnum][(2*i)*257+(2*j-1)];
                if (i>0)
                    fs2 >> table[partnum][(2*i-1)*257+(2*j-1)];
            }
        }

    fs1.close();
    fs2.close();
    fs3.close();*/

    //generate nonminimal cos and sin matrices
    std::cerr << "OK fillStates ";
    std::vector<double> chmat(dim*dim);
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if ((H.GetState(i).ParticleNumber()%2)!=(H.GetState(j).ParticleNumber()%2))
            {
                chmat[i*dim+j] = 0.0;
            } else
                chmat[i*dim+j]=zmu[l-1]*fmass*l*CosMxEl(H.GetState(i),H.GetState(j),table,maxnum);
        }
    }
    std::cerr << "OK AssembleMx: ";

    // generate nonminimal diagonal matrix H0
  /*  std::vector<double> energymat(dim*dim);*/
    for(int i=0; i< dim; i++)
        chmat[i*(dim+1)]+=H.Energy0(i);

    // set up full Hamiltonian!

    //diagonalize
    std::cerr << "OK Diagonalize:\n ";

    double* mx = (double*)calloc(chmat.size(),sizeof(double));
    double* evals = (double*)calloc(dim,sizeof(double));

 //LAPACK diagonalization

    std::copy(chmat.begin(),chmat.end(),mx);
    int info;
    char cV= 'N';
    char cU= 'L';
    double worksize;
    int lwork=-1;
    LAPACK_dsyev(&cV,&cU,&dim,mx,&dim,evals,&worksize,&lwork,&info);
    lwork=worksize;

    double* work=(double*)calloc(lwork,sizeof(double));

    LAPACK_dsyev(&cV,&cU,&dim,mx,&dim,evals,work,&lwork,&info);

    free(work);
    free(mx);
//    std::cout << "Last eval=" << evals[0] << '\n';

 //PRIMME diagonalization

 /*   sparse_matrix mx_s;
    dense2sparse(dim,chmat,mx_s);
    chmat.clear();
    mx_s.sort();
    double* evecs =(double*)calloc(dim*10,sizeof(double));
    comp_evals(dim,0.1, 7, mx_s, evals, evecs);
    free(evecs);
*/
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

    free(evals);

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
}



int main(int argc, char* argv[])
{
    if (argc != 5)
    {
        std::cout << "usage: " << argv[0] << "<directory_in> <partnum_in> <parity_in> <mass_sign>\n";
        return 0;
    }
 /*   std::fstream mattest;
    mattest.open("mattest",std::ios::out | std::ios::trunc);
    for (int i=0; i < 25; i++)
    {
        for (int j =0; j < 25; j++)
                mattest << MinCosMat(i,j,20.0,24) << '\t';
        mattest << '\n';
    }
    mattest.close();*/

    double M = 1.0;

    short PBC = 1;
    short partnum = atoi(argv[2]);
    partnum_in = partnum;
    double E0max = partnum*2.0;
    parity_in = atoi(argv[3]);
    char* dir =argv[1];
    double ms = atoi(argv[4]);
    std::vector<double> fmassvec {ms*0.01,ms*0.02,ms*0.04,ms*0.06,ms*0.07,ms*0.08,ms*0.09,ms*0.1,ms*0.2,ms*0.4,ms*0.7,ms*1.0};
    std::vector<double> thetavec {ms*pi/12.0,ms*pi/6.0,ms*pi/4.0,ms*pi/3.0,ms*5.0/12.0*pi,ms*pi/2.0,ms*7.0/12.0*pi,ms*2.0/3.0*pi,ms*3.0/4.0*pi,ms*5.0/6.0*pi,ms*11.0/12.0*pi,ms*pi};

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
    //HilbertSpace He;
    //HilbertSpace Ho;

    H.setParams(1.0*pi,1.0,PBC);
    //He.setParams(2.0*pi,1.0,PBC);
    //Ho.setParams(2.0*pi,1.0,PBC);
    HilbertState vacuum;

    E0max = basis(H,vacuum,0, E0max+0.1,partvec[partnum],parity_in,0);
    //E0max = basis(He,vacuum,0, E0max+0.1,partvec[partnum],0,0); //partvec[partnum]

    std::fstream enfile;
    enfile.open("energies",std::ios::out | std::ios::app);
    enfile << dir << " " << partnum << " " << parity_in << " " << ms << " " << E0max << ' ';
 //   E0max = basis(Ho,vacuum,0, partnum*2.0+0.1,partvec[partnum],1,0);
 //   enfile << E0max << '\n';
    enfile.close();

    if (E0max<0.0)
        return -1;
    std::cout << H.NumStates() << " states generated\n";
 //   std::cout << He.NumStates() << '+' << Ho.NumStates() << " states generated\n";
 //   for (int t=0; t <= 11; t++)
 //   {
        for (int i=8; i <= 8; i++)
        {
            std::string fname = std::string(dir)+"/m"+std::to_string(i+1)+"/thmL";
            //std::string fname = std::string(dir)+"/th"+std::to_string(t+1)+"/m"+std::to_string(i+1)+"/thmL";
            for (int l=21; l <= 40; l++)
            {
                H.setParams(l,M,PBC);
              //  E0max = basis(H,vacuum,0, 10,partvec[partnum],parity_in,0);
             //   He.setParams(l,M,PBC);
             //   Ho.setParams(l,M,PBC);
                calculateSpectrum(H, fmassvec[i], partnum, l, fname);
 //               calculateSpectrumTh(He,Ho, fmassvec[i], thetavec[t], partnum, l, fname);
            }
        }
//    }

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
