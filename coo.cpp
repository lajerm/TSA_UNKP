#include "coo.h"
#include <iostream>

#include "phi4.h"

bool COO_comp(COO c1, COO c2)
{
    if (c1.row_ind < c2.row_ind)
    {
        return true;
    }
    if (c1.row_ind == c2.row_ind)
    {
        if (c1.col_ind < c2.col_ind)
            return true;
    }
    return false;
}

sparse_matrix::sparse_matrix( std::vector<COO>& mx_in)
{
    matrix = mx_in;
    std::sort(matrix.begin(),matrix.end(),COO_comp);
}

sparse_matrix::sparse_matrix(const std::vector<sparse_matrix>& spvec_in)
{
    matrix.clear();
    for (int i = 0; i < spvec_in.size(); i++)
    {
        matrix.insert(matrix.end(),spvec_in[i].matrix.begin(),spvec_in[i].matrix.end());
    }
    std::sort(matrix.begin(),matrix.end(),COO_comp);
    int i = 0;
    if (matrix.size() > 0)
    {
    double val = matrix[0].val;
    COO coo;

    for (int j = 1; j < matrix.size(); j++)
    {
        if ((matrix[j].col_ind == matrix[j-1].col_ind)&&(matrix[j].row_ind==matrix[j-1].row_ind))
        {
            val += matrix[j].val;
        } else
        {
            coo.row_ind = matrix[j-1].row_ind;
            coo.col_ind = matrix[j-1].col_ind;
            coo.val = val;
            if (val != 0.0)
            {
                matrix[i] = coo;
                i++;
            }
            val = matrix[j].val;
        }
    }

    if (val != 0.0)
    {
        coo.row_ind = matrix[matrix.size()-1].row_ind;
        coo.col_ind = matrix[matrix.size()-1].col_ind;
        coo.val = val;
        matrix[i] = coo;
        i++;
    }
        matrix.resize(i);
    }
}
double sparse_matrix::add(int i, int j, double num)
{
    COO coo;
    if (i==j)   //diagonal elements
    {
        coo.row_ind = i;
        coo.col_ind = j;
        coo.val = num;
        matrix.push_back(coo);
        return sqr(num);
    }
    else    //off-diagonal elements
    if (num != 0.0)
    {
        coo.row_ind = i;
        coo.col_ind = j;
        coo.val = num;
        matrix.push_back(coo);
        coo.row_ind = j;
        coo.col_ind = i;
        matrix.push_back(coo);
        return 2.0*sqr(num);
    }
    return 0.0;
}

void sparse_matrix::add_disjunct_part(int i,int j, const sparse_matrix &spvec) //i: row shift, j: column shift
{
    COO coo;
    for (int k = 0; k < spvec.size(); k++) //this must be refined: even/odd submatrix sizes vary with i,j
    {
        coo.row_ind = spvec.matrix[k].row_ind+i;
        coo.col_ind = spvec.matrix[k].col_ind+j;
        coo.val = spvec.matrix[k].val;

        matrix.push_back(coo);
        if (i!=j)
        {
            coo.row_ind = spvec.matrix[k].col_ind+j;
            coo.col_ind = spvec.matrix[k].row_ind+i;
            matrix.push_back(coo);
        }
    }
}

sparse_matrix sparse_matrix::mul (double c)
{
    sparse_matrix temp = *this;
    for (int i = 0; i < matrix.size(); i++)
        temp.matrix[i].val *= c;
        return temp;
}


double sparse_matrix::COO_trunc(int index)
{
    double fnorm = 0;
    int i = 0;
    for (int j = 0; j < matrix.size(); j++)
    {
        if ((matrix[j].col_ind < index)&&(matrix[j].row_ind < index))
        {
            matrix[i] = matrix[j];
            fnorm += sqr(matrix[i].val);
            i++;
        }
    }
    matrix.resize(i);
    return fnorm;
}
