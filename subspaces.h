#ifndef __SUBSPACES_H__
#define __SUBSPACES_H__

#include <vector>
#include <functional>
#include "coo.h"

int directProdDD(int k1, int l1, int k2, int l2, const std::vector<double> &d1_in, const std::vector<double> &d2_in, std::vector<double> &d_out);
int directProdDS(int k1, int l1, int k2, int l2, const std::vector<double> &d1_in, const sparse_matrix &s2_in, sparse_matrix &s_out);
int directProdSS(int k2, int l2, const sparse_matrix &s1_in, const sparse_matrix &s2_in, sparse_matrix &s_out);

class basic_operator_dense
{
    int sector1;
    int sector2;
    std::vector<double> data;
public:
    basic_operator_dense(int s1, int s2, const std::vector<double> &d_in) : sector1(s1), sector2(s2), data(d_in) {}
    int getSector1() {return sector1;}
    int getSector2() {return sector2;}
    int get(std::vector<double> &data_out) { data_out = data; return 0; }
};

class basic_operator_sparse
{
    int sector1;
    int sector2;
    sparse_matrix data;
public:
    basic_operator_sparse(int s1, int s2, const sparse_matrix & s_in) : sector1(s1), sector2(s2), data(s_in) {}

    int getSector1() {return sector1;}
    int getSector2() {return sector2;}
    int get(sparse_matrix &data_out) {data_out = data; return 0; }
};

typedef std::vector<basic_operator_dense> operator_dense;
typedef std::vector<basic_operator_sparse> operator_sparse;

template<class T> class subSpace
{
    std::vector<T> operators;
    std::vector<std::vector<int> > sectorInfo;
    std::vector<int> sectorDims;
public:
    subSpace(const std::vector<T> &op_in, const std::vector<std::vector<int> > secInf_in, const std::vector<int> secDims_in) : operators(op_in), sectorInfo(secInf_in), sectorDims(secDims_in) {}
    int getSectorInfo(int i, std::vector<int> &info) {info = sectorInfo[i]; return 0; }
    short numSectors() {return sectorInfo.size(); }
    int getDim(int i) { return sectorDims[i]; }
    int getOperator(int i, T &op_out) {op_out = operators[i]; return 0; }
};

class StateSpace
{
    std::vector<subSpace<operator_dense> > denseSubspaces;
    std::vector<subSpace<operator_sparse> > sparseSubspaces;
    std::vector<std::vector<int> > allowedSubspaces;

    std::vector<std::vector<std::vector<int> > > allowedStates;

public:
    StateSpace(const std::vector<subSpace<operator_dense> > &dsS_in, const std::vector<subSpace<operator_sparse> > &ssS_in) : denseSubspaces(dsS_in), sparseSubspaces(ssS_in) {}
    int generateAllowedSubspaces(bool (*compare)(const std::vector<std::vector<int> >&));
    int nonzeroMatrix(int subspaceID, int operatorID, int sectorID1, int sectorID2, short &transpose);
    int DirectProductDD(int subspaceID, const std::vector<double> &vec_in, int operatorID, int sectorID1, int sectorID2, std::vector<double> &vec_out, int &dimR, int &dimC);
    int DirectProductDS(const std::vector<double> &vec_in, int operatorID, int sectorID1, int sectorID2, sparse_matrix &sparse_out, int &dimR, int &dimC);
    int DirectProductSS(int subspaceID, const sparse_matrix &sparse_in, int operatorID, int sectorID1, int sectorID2, sparse_matrix &sparse_out, int &dimR, int &dimC);

    int numAllowedSubspaces() {return allowedSubspaces.size();}
    int numDense() {return denseSubspaces.size();}
    int numSparse() {return sparseSubspaces.size();}
    int getAllowedSubspaceDim(int i)
    {
        int dim=1;
        for (int j=0; j < denseSubspaces.size(); j++) dim *= denseSubspaces[j].getDim(allowedSubspaces[i][j]);
        for (int j=0; j < sparseSubspaces.size(); j++) dim *=sparseSubspaces[j].getDim(allowedSubspaces[i][j+denseSubspaces.size()]);
        return dim;
    }

    int generateAllowedStates(std::function<int(const std::vector<int> &, std::vector<std::vector<int> > &)> genStates )
    {
        for (int i=0; i < allowedSubspaces.size(); i++)
        {
            std::vector<std::vector<int> > temp;
            genStates(allowedSubspaces[i],temp);
            allowedStates.push_back(temp);
            if (temp.size()>0) temp.clear();
        }
    }
    int getAllowedStates(std::vector<std::vector<std::vector<int> > > &states_out) {states_out = allowedStates; return 0;}
    int getKeptAllowedSPDim(int i) {return allowedStates[i].size();}

    std::vector<double> getDenseOperator(int subspaceID, int operatorID, int sectorID1, int sectorID2, int &dim1, int &dim2);

};

class HamiltonianDense
{
    std::vector<double> coefficients;
    std::vector<std::vector<int> > interaction;
    StateSpace *StateSp;
    std::vector<double> matrix;
public:
    HamiltonianDense(const std::vector<double> &coefs_in, const std::vector<std::vector<int> > &interact_in, StateSpace *StateSp_in): coefficients(coefs_in), interaction(interact_in), StateSp(StateSp_in) {}
    int generateHamiltonian();

    int extractMatrix(std::vector<double> &mx_out, int &N) {mx_out = matrix;  N=0; for (int i=0; i<StateSp->numAllowedSubspaces(); i++) N+= StateSp->getKeptAllowedSPDim(i); return 0;}
    int clear() {coefficients.clear(); interaction.clear(); matrix.clear();}
};

class HamiltonianSparse
{
    std::vector<double> coefficients;
    std::vector<std::vector<int> > interaction;
    StateSpace *StateSp;
    sparse_matrix matrix;
public:
    HamiltonianSparse(const std::vector<double> coefs_in, const std::vector<std::vector<int> > &interact_in, StateSpace *StateSp_in): coefficients(coefs_in), interaction(interact_in), StateSp(StateSp_in) {}
    int generateHamiltonian();

    int extractMatrix(sparse_matrix &mx_out) {mx_out = matrix; return 0;}

};


#endif // __SUBSPACES_H__
