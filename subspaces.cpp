
#include "subspaces.h"
#include <iostream>
//testing purposes

int StateSpace::generateAllowedSubspaces(bool (*compare)(const std::vector<std::vector<int> >&))
{
    if (allowedSubspaces.size()>0) allowedSubspaces.clear();

    std::vector<int> markers;
    int nD = denseSubspaces.size();
    int nS = sparseSubspaces.size();
    std::vector<int> maxSector;
    std::vector<std::vector<int> > sectorInfoVec;

    markers.resize(nD+nS,0);

    for (int i = 0; i < nD; i++)
    {
        std::vector<int> tempSectorInfo;
        maxSector.push_back(denseSubspaces[i].numSectors());
        denseSubspaces[i].getSectorInfo(markers[i],tempSectorInfo);
        sectorInfoVec.push_back(tempSectorInfo);
        tempSectorInfo.clear();
    }
    for (int i = 0; i < nS; i++)
    {
        std::vector<int> tempSectorInfo;
        maxSector.push_back(sparseSubspaces[i].numSectors());
        sparseSubspaces[i].getSectorInfo(markers[i+nD],tempSectorInfo);
        sectorInfoVec.push_back(tempSectorInfo);
        tempSectorInfo.clear();
    }
    if (compare(sectorInfoVec))
        allowedSubspaces.push_back(markers);

    sectorInfoVec.clear();

    while (markers[nD+nS-1] < maxSector[nD+nS-1])
    {
        markers[0]++;
        int i = 0;
        while ((markers[i]==maxSector[i])&&(i<(nD+nS-1)))
        {
            i++;
            if (i < (nD+nS))
                markers[i]++;
            for (int j = 1; j <= i; j++)
                markers[i-j] = 0;
        }
        if (markers[nD+nS-1] < maxSector[nD+nS-1])
        {
            for (int i = 0; i < nD; i++)
            {
                std::vector<int> tempSectorInfo;
                denseSubspaces[i].getSectorInfo(markers[i],tempSectorInfo);
                sectorInfoVec.push_back(tempSectorInfo);
                tempSectorInfo.clear();
            }
            for (int i = 0; i < nS; i++)
            {
                std::vector<int> tempSectorInfo;
                sparseSubspaces[i].getSectorInfo(markers[i+nD],tempSectorInfo);
                sectorInfoVec.push_back(tempSectorInfo);
                tempSectorInfo.clear();
            }
            if (compare(sectorInfoVec))
                allowedSubspaces.push_back(markers);

            sectorInfoVec.clear();
        }
    }
/*    std::cout << "ALLOWED SUBSPACES TEST\n";
    for (int i = 0; i < allowedSubspaces.size(); i++)
    {
        for (int j =0; j < allowedSubspaces[i].size(); j++) std::cout << allowedSubspaces[i][j] << ' ';
        std::cout << '\n';
    }*/

    return 0;
}

int StateSpace::nonzeroMatrix(int subspaceID, int operatorID, int sectorID1, int sectorID2, short &transpose)
{
    int index = -1;
    if (subspaceID < denseSubspaces.size())     // we are in the dense part of the direct product
    {
        operator_dense Ov;
        denseSubspaces[subspaceID].getOperator(operatorID,Ov);
        for (int m =0; m < Ov.size(); m++) //check if the given operator has a nonzero matrix in the given sector
        {
            if ((Ov[m].getSector1()==allowedSubspaces[sectorID1][subspaceID])&&
                    (Ov[m].getSector2()==allowedSubspaces[sectorID2][subspaceID]))
                    {
                        index=m;
                        transpose = 0;
                    }
            if ((Ov[m].getSector1()==allowedSubspaces[sectorID2][subspaceID])&&
                    (Ov[m].getSector2()==allowedSubspaces[sectorID1][subspaceID]))
                    {
                        index=m;
                        transpose = 1;
                    }
        }
        Ov.clear();
    }
    else // we are in the sparse part of the direct product
    {
        operator_sparse Ov;
        sparseSubspaces[subspaceID-denseSubspaces.size()].getOperator(operatorID,Ov);
        for (int m =0; m < Ov.size(); m++) //check if the given operator has a nonzero matrix in the given sector
        {
            if ((Ov[m].getSector1()==allowedSubspaces[sectorID1][subspaceID])&&
                    (Ov[m].getSector2()==allowedSubspaces[sectorID2][subspaceID]))
                    {
                        index=m;
                        transpose = 0;
                    }
            if ((Ov[m].getSector1()==allowedSubspaces[sectorID2][subspaceID])&&
                    (Ov[m].getSector2()==allowedSubspaces[sectorID1][subspaceID]))
                    {
                        index=m;
                        transpose = 1;
                    }
        }
    }
    return index;
}

int directProdDD(int k1, int l1, int k2, int l2, const std::vector<double> &d1_in, const std::vector<double> &d2_in, std::vector<double> &d_out, short transpose)
{
    for (int i1 = 0; i1 < k1; i1++)
    {
        for (int i2 = 0; i2 < k2; i2++ )
        {
            for (int j1 = 0; j1 < l1; j1++)
            {
                for (int j2 = 0; j2 < l2; j2++)
                {
                    if (transpose == 0)
                    {
                        d_out.push_back(d1_in[i1*l1+j1]*d2_in[i2*l2+j2]);
                    } else
                    {
                        d_out.push_back(d1_in[i1*l1+j1]*d2_in[j2*l2+i2]);
                    }

                }
            }
        }
    }
    return 0;
}

int directProdDS(int k1, int l1, int k2, int l2, const std::vector<double> &d1_in, const sparse_matrix &s2_in, sparse_matrix &s_out, short transpose)
{
    COO coo;
    for (int i = 0; i < d1_in.size(); i++)
    {
        for (int j = 0; j < s2_in.size(); j++)
        {
            coo=s2_in.get(j);
            if (transpose==1)
            {
                int t = coo.col_ind;
                coo.col_ind = coo.row_ind;
                coo.row_ind = t;
            }
            coo.row_ind += (i/k1)*k2;
            coo.col_ind += (i%k1)*l2;
            coo.val *= d1_in[i];
            if (coo.val != 0.0)
                s_out.PushBack(coo);
        }
        s_out.sort();
    }
    return 0;
}

int directProdSS(int k2, int l2, const sparse_matrix &s1_in, const sparse_matrix &s2_in, sparse_matrix &s_out, short transpose)
{
    COO coo1, coo2;
    for (int i = 0; i < s1_in.size(); i++)
    {
        coo1 = s1_in.get(i);

        for (int j = 0; j < s2_in.size(); j++)
        {
            coo2 = s2_in.get(j);
            if (transpose==1)
            {
                int t = coo2.col_ind;
                coo2.col_ind = coo2.row_ind;
                coo2.row_ind = t;
            }
            coo2.row_ind += coo1.row_ind*k2;
            coo2.col_ind += coo1.col_ind*l2;
            coo2.val *=coo1.val;
            if (coo2.val != 0.0)
                s_out.PushBack(coo2);
        }
    }
    s_out.sort();
    return 0;
}

int StateSpace::DirectProductDD(int subspaceID, const std::vector<double> &vec_in, int operatorID, int sectorID1, int sectorID2, std::vector<double> &vec_out, int &dimR, int &dimC)
// l, denseTemp, interaction[k][l],i,j,denseTemp2, dimR, dimC
{
    short transpose = 0;
    int index=nonzeroMatrix(subspaceID,operatorID,sectorID1,sectorID2,transpose);
    operator_dense Ov;
    denseSubspaces[subspaceID].getOperator(operatorID,Ov);
    std::vector<double> op;
    Ov[index].get(op);
    directProdDD(dimR,dimC,denseSubspaces[subspaceID].getDim(allowedSubspaces[sectorID1][subspaceID]),denseSubspaces[subspaceID].getDim(allowedSubspaces[sectorID2][subspaceID]),vec_in,
                op,vec_out,transpose);
    if (!transpose)
    {
        dimR *= denseSubspaces[subspaceID].getDim(allowedSubspaces[sectorID1][subspaceID]);
        dimC *= denseSubspaces[subspaceID].getDim(allowedSubspaces[sectorID2][subspaceID]);
    } else
    {
        dimR *= denseSubspaces[subspaceID].getDim(allowedSubspaces[sectorID2][subspaceID]);
        dimC *= denseSubspaces[subspaceID].getDim(allowedSubspaces[sectorID1][subspaceID]);

    }
        op.clear();
    return 0;
}

std::vector<double> StateSpace::getDenseOperator(int subspaceID, int operatorID, int sectorID1, int sectorID2, int &dim1, int &dim2)
{
    short transpose = 0;
    int index=nonzeroMatrix(subspaceID,operatorID,sectorID1,sectorID2,transpose);
    operator_dense Ov;
    denseSubspaces[subspaceID].getOperator(operatorID,Ov);
    std::vector<double> op;
    Ov[index].get(op);
    dim1 = denseSubspaces[subspaceID].getDim(allowedSubspaces[sectorID1][subspaceID]);
    dim2 = denseSubspaces[subspaceID].getDim(allowedSubspaces[sectorID2][subspaceID]);
    if (transpose)
    {
        std::vector<double> op_(dim1*dim2);
        for (int i=0; i < dim1; i++)
        {
            for (int j=0; j < dim2; j++)
                op_[j*dim1+i]=op[i*dim2+j];
        }
        op=op_;
        dim2=dim1;
        dim1=op_.size()/dim1;
    }
    return op;
}

int StateSpace::DirectProductDS(const std::vector<double> &vec_in, int operatorID, int sectorID1, int sectorID2, sparse_matrix &sparse_out, int &dimR, int &dimC)
// denseTemp, interaction[k][StateSp->denseSubspaces.size()], i,j, sparse_temp, dimR,dimC
{
    short transpose = 0;
    int index=nonzeroMatrix(numDense(),operatorID,sectorID1,sectorID2,transpose);
    operator_sparse Ov;
    sparseSubspaces[0].getOperator(operatorID,Ov);
    sparse_matrix op;
    Ov[index].get(op);
    directProdDS(dimR,dimC,sparseSubspaces[0].getDim(allowedSubspaces[sectorID1][numDense()]),
                sparseSubspaces[0].getDim(allowedSubspaces[sectorID2][numDense()]),vec_in,
                op,sparse_out,transpose);
    if (!transpose)
    {
        dimR *= sparseSubspaces[0].getDim(allowedSubspaces[sectorID1][numDense()]);
        dimC *= sparseSubspaces[0].getDim(allowedSubspaces[sectorID2][numDense()]);
    } else
    {
        dimR *= sparseSubspaces[0].getDim(allowedSubspaces[sectorID2][numDense()]);
        dimC *= sparseSubspaces[0].getDim(allowedSubspaces[sectorID1][numDense()]);
    }

    op.clear();
}
int StateSpace::DirectProductSS(int subspaceID, const sparse_matrix &sparse_in, int operatorID, int sectorID1, int sectorID2, sparse_matrix &sparse_out, int &dimR, int &dimC)
// l+denseSubspaces.size(), sparse_temp, interaction[k][l+StateSp->denseSubspaces.size()], i,j, sparse_temp2, dimR,dimC
{
    short transpose = 0;
    int index=nonzeroMatrix(subspaceID,operatorID,sectorID1,sectorID2,transpose);
    operator_sparse Ov;
    sparseSubspaces[subspaceID-numDense()].getOperator(operatorID,Ov);
    sparse_matrix op;
    Ov[index].get(op);
    directProdSS(sparseSubspaces[subspaceID-numDense()].getDim(allowedSubspaces[sectorID1][subspaceID]),
    sparseSubspaces[subspaceID-numDense()].getDim(allowedSubspaces[sectorID2][subspaceID]),sparse_in,
    op,sparse_out,transpose);
    if (!transpose)
    {
        dimR *= sparseSubspaces[subspaceID-numDense()].getDim(allowedSubspaces[sectorID1][subspaceID]);
        dimC *= sparseSubspaces[subspaceID-numDense()].getDim(allowedSubspaces[sectorID2][subspaceID]);
    } else
    {
        dimR *= sparseSubspaces[subspaceID-numDense()].getDim(allowedSubspaces[sectorID2][subspaceID]);
        dimC *= sparseSubspaces[subspaceID-numDense()].getDim(allowedSubspaces[sectorID1][subspaceID]);
    }

    op.clear();
}

 int HamiltonianSparse::generateHamiltonian( /* függvényargumentum feltétel az állapotokra */)
{
    int rowShift = 0;
    for (int i = 0; i < StateSp->numAllowedSubspaces(); i++)     // rows of full Hamiltonian hypermatrix
    {

        int colShift = rowShift;
        for (int j = i; j < StateSp->numAllowedSubspaces(); j++)     // columns of full Hamiltonian hypermatrix
        {
            std::vector<sparse_matrix> sparsePart;       // store sparse part of direct product operator sum
            // make a list of relevant interactions
            int dimR = 1;
            int dimC = 1;
            if (sparsePart.size()>0) sparsePart.clear();
            for (int k = 0; k < interaction.size(); k++)    // go through interaction terms
            {
                bool flag = true;   // this flag will turn false if the opertator is zero in the given sector
                int l = 0;
                while ((l<interaction[k].size())&&flag) // go through the factors of an interaction term
                {
                    short transpose = 0;
                    if (StateSp->nonzeroMatrix(l,interaction[k][l],i,j,transpose)<0) flag = false;
                    l++;
                }
                if (flag) // the operator of interaction[k] does have nonzero elements in the sector pair i,j
                {
                    dimR = 1;             //in fact we could construct the direct product matrix here
                    dimC = 1;
                    std::vector<double> denseTemp{1.0};
                    std::vector<double> denseTemp2;
                    for (int l = 0; l < StateSp->numDense(); l++)   // factor up the dense part
                    {
                        StateSp->DirectProductDD(l, denseTemp, interaction[k][l],i,j,denseTemp2,dimR,dimC);
                        denseTemp=denseTemp2;
                        denseTemp2.clear();
                    }
                    sparse_matrix sparse_temp;
                    StateSp->DirectProductDS(denseTemp, interaction[k][StateSp->numDense()], i,j, sparse_temp,dimR,dimC);  // factor the dense result with first sparse term
                    sparse_matrix sparse_temp2;
                    for (int l = 1; l < StateSp->numSparse(); l++)   // factor up the sparse part
                    {
                        StateSp->DirectProductSS(l+StateSp->numDense(), sparse_temp, interaction[k][l+StateSp->numDense()], i,j, sparse_temp2,dimR,dimC);
                        sparse_temp = sparse_temp2;
                        sparse_temp2.clear();
                    }

                    sparsePart.push_back(sparse_temp.mul(coefficients[k]));  // collect all operators in the sector
                }
            }
            // dimR .és dimC helyére a dimenziók szorzatai kellenek!
            //std::cout << dimR << ',' << dimC << '\t';
            std::cout << StateSp->getAllowedSubspaceDim(i) << ',' << StateSp->getAllowedSubspaceDim(j) << '\n';
            matrix.add_disjunct_part(rowShift,colShift,sparse_matrix(sparsePart));  // for offdiagonal submatrices, this adds the transpose, too

            colShift += StateSp->getAllowedSubspaceDim(j);
        }
        rowShift += StateSp->getAllowedSubspaceDim(i);
    }
    return 0;
}

int HamiltonianDense::generateHamiltonian()
{
    int rowShift = 0;
    std::vector<std::vector<std::vector<int> > > stateInfo;
    StateSp->getAllowedStates(stateInfo);
    int N =0;

    for (int i=0; i<stateInfo.size(); i++) N += stateInfo[i].size();
    matrix.resize(N*N);

    for (int i = 0; i < StateSp->numAllowedSubspaces(); i++)     // rows of full Hamiltonian hypermatrix
    {

        int colShift = rowShift;
        for (int j = i; j < StateSp->numAllowedSubspaces(); j++)     // columns of full Hamiltonian hypermatrix
        {
            // make a list of relevant interactions
            for (int k = 0; k < interaction.size(); k++)    // go through interaction terms
            {
                bool flag = true;   // this flag will turn false if the opertator is zero in the given sector
                int l = 0;
                while ((l<interaction[k].size())&&flag) // go through the factors of an interaction term
                {
                    short transpose = 0;
                    if (StateSp->nonzeroMatrix(l,interaction[k][l],i,j,transpose)<0) flag = false;
                    l++;
                }
                if (flag) // the operator of interaction[k] does have nonzero elements in the sector pair i,j
                {
                    std::vector<std::vector<double> > Ov(StateSp->numDense());
                    std::vector<int> widths(StateSp->numDense());
                    int dim1;
                    for (int l = 0; l < StateSp->numDense(); l++)
                    {
                         Ov[l]=StateSp->getDenseOperator(l,interaction[k][l],i,j,dim1,widths[l]);
                    }


                    for (int l=0; l < stateInfo[i].size(); l++)
                        for (int m=0; m < stateInfo[j].size(); m++)
                        {
                            double MX = 1.0;

                            for (int n=0; n < StateSp->numDense(); n++)
                            {
                               MX *= Ov[n][stateInfo[i][l][n]*widths[n]+stateInfo[j][m][n]];

                            }

                            matrix[(rowShift+l)*N+colShift+m] += coefficients[k]*MX;
                            if (rowShift != colShift)
                                matrix[(colShift+m)*N+rowShift+l] += coefficients[k]*MX;
                        }
                }
            }
            colShift += stateInfo[j].size();
        }
        rowShift += stateInfo[i].size();
    }
    return 0;
}
