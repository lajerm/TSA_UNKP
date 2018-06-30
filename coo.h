#ifndef __COO_H__
#define __COO_H__

#include <vector>
#include <algorithm>
// #include "utils.h"


struct COO
{
    double val;
    int col_ind;
    int row_ind;
};

bool COO_comp(COO c1, COO c2);

class sparse_matrix
{
    std::vector<COO> matrix;
public:
    sparse_matrix() { }
    sparse_matrix( std::vector<COO>& mx_in);
    sparse_matrix(const std::vector<sparse_matrix>& spvec_in);

    ~sparse_matrix() {matrix.clear();}

    void obtainMx(std::vector<COO> & mx_out) const {mx_out = matrix;}

    COO get(int i) const  {return matrix[i];}
    int size() const {return matrix.size();}
    void set(int i, const COO& coo_in) {matrix[i] = coo_in;}
    double add(int i, int j, double num);
    void PushBack(COO coo) {matrix.push_back(coo); }
    void clear() {matrix.clear();}
    void sort() {std::sort(matrix.begin(), matrix.end(), COO_comp); }

    void add_disjunct_part(int i,int j,const sparse_matrix &spvec);
    sparse_matrix mul (double c);

    double COO_trunc(int index);

};


#endif // __COO_H__
