#ifndef _CSC_HPP_
#define _CSC_HPP_

#include <vector>
#include <iostream>
#include <stdlib.h>

#include "sparse.hpp"

using namespace std;

template< typename ValueType, typename IndicesType, int base>
void CooToCsc(CooMatrix<ValueType,IndicesType,base> &R, CscMatrix<ValueType,IndicesType,base> &L);

//! Compressed row storage format sparse matrix
/*!
 * @tparam ValueType the type of the entries of the matrix, for the moment only double is supported
 * @tparam IndicesType the type of vector that stores row/column indices, default is size_t
 * @tparam base the base of the vector indices, 0 for c-based, 1 for fortran-based
 */
template<
    typename ValueType = double,
    typename IndicesType = int,
    int base = 0
    >
class CscMatrix : public SparseMatrix<ValueType,IndicesType,base> {
    protected:
        //! Creates an compressed row storage sparse matrix with vectors as reference
        /*!
         * @param m number of rows
         * @param n number of columns
         * @param nz number of nonzero entries
         * @param iptr the row pointer vector
         * @param jcn the column indices vector
         * @param val the entries vector
         * @param reference the vectors in the matrix are the real ones 
         * @param switchToBase if @base == 0 then apply a -1 to indices, if @base == 1 then apply a +1
         */
        void init( IndicesType m,  IndicesType n,  IndicesType nz,
                 IndicesType* irn,  IndicesType* jptr,  ValueType* val,
                 bool reference, bool swithToBase)
        {
            this->m_ = m;
            this->n_ = n;
            this->nz_ = nz;
            this->isReference = reference;
            if (this->isReference) {
                this->val_  = val;
                this->rows_ = irn;
                this->cols_ = jptr;
            } else {
                this->val_ = new ValueType[this->nz_];
                this->rows_ = new IndicesType[this->nz_];
                this->cols_ = new IndicesType[this->n_ + 1];
                if (swithToBase) {
                    ValueType b;
                    if (base == 0) {
                        b = -1;
                    } else {
                        b = 1;
                    }
                    for (IndicesType t = 0; t < this->nz_; t++) {
                        this->val_[t] = val[t];
                        this->rows_[t] = irn[t] + b;
                    }
                    for (IndicesType t = 0; t < this->m_; t++) {
                        this->cols_[t] = jptr[t] + b;
                    }
                } else {
                    std::copy(val, val + this->nz_, this->val_);
                    std::copy(irn, irn + this->nz_, this->rows_);
                    std::copy(jptr, jptr + this->n_ + 1, this->cols_);
                }
            }
            this->isSymmetric = false;
        }
    public:
        /*!
         * Basic constructor, initialize everything to zero
         */
        CscMatrix()
        {
            this->m_ = 0;
            this->n_ = 0;
            this->nz_ = 0;
            this->val_ = 0;
            this->rows_ = 0;
            this->cols_ = 0;
            this->isReference = false;
            this->isSymmetric = false;
        }

        //! Creates an empty compressed row storage sparse matrix
        /*!
         * @param m number of rows
         * @param n number of columns
         * @param nz number of nonzero entries
         */
        CscMatrix( IndicesType m,  IndicesType n,  IndicesType nz)
        {
            this->m_ = m;
            this->n_ = n;
            this->nz_ = nz;
            this->val_ = new ValueType[this->nz_];
            this->rows_ = new IndicesType[this->nz_];
            this->cols_ = new IndicesType[this->n_ + 1];
            this->isReference = false;
            this->isSymmetric = false;
            this->cols_[this->n_] = this->nz_;
        }


        //! Creates an compressed row storage sparse matrix
        /*!
         * @param m number of rows
         * @param n number of columns
         * @param nz number of nonzero entries
         * @param iptr the row pointer vector
         * @param jcn the column indices vector
         * @param val the entries vector
         */
        CscMatrix( IndicesType m,  IndicesType n,  IndicesType nz,
                 IndicesType* irn,  IndicesType* jptr,  ValueType* val)
        {
            this->init(m, n, nz, irn, jptr, val, false, false);
        }

        //! Creates an compressed row storage sparse matrix
        /*!
         * @param m number of rows
         * @param n number of columns
         * @param nz number of nonzero entries
         * @param iptr the row pointer vector
         * @param jcn the column indices vector
         * @param val the entries vector
         * @param reference the vectors in the matrix are the real ones 
         */
        CscMatrix( IndicesType m,  IndicesType n,  IndicesType nz,
                 IndicesType* irn,  IndicesType* jptr,  ValueType* val,
                 bool reference)
        {
            this->init(m, n, nz, irn, jptr, val, reference, false);
        }

        //! Creates an compressed row storage sparse matrix
        /*!
         * @param m number of rows
         * @param n number of columns
         * @param nz number of nonzero entries
         * @param iptr the row indices vector
         * @param jcn the column indices vector
         * @param val the entries vector
         * @param reference the vectors in the matrix are the real ones 
         * @param switchToBase if @base == 0 then apply a -1 to indices, if @base == 1 then apply a +1
         */
        CscMatrix( IndicesType m,  IndicesType n,  IndicesType nz,
                 IndicesType* irn,  IndicesType* jptr,  ValueType* val,
                 bool reference, bool swithToBase)
        {
            this->init( m,  n,  nz, irn, jptr, val, reference, swithToBase);
        }

        //! Creates a compressed row storage sparse matrix from a coordinate storage matrix
        /*!
         * @param C the coordinate storage matrix to be converted
         */
        CscMatrix(CooMatrix<ValueType,IndicesType,base> &C)
        {
            this->m_ = C.nRows();
            this->n_ = C.nCols();
            this->nz_ = C.nnz();
            this->val_ = new ValueType[this->nz_];
            this->rows_ = new IndicesType[this->nz_];
            this->cols_ = new IndicesType[this->n_ + 1];
            this->isReference = false;
            this->isSymmetric = false;
            this->cols_[this->n_] = this->nz_;

            CooToCsc(C, *this);
        }

        //! Assign the right-hand side matrix into the left-hand side (reference)
        /*!
         * Returns a reference to the left-hand side matrix, if you want a copy
         * use the constructor. If the matrix contained data, it will be deleted.
         *
         * @return a reference to the left-hand side
         */
        CscMatrix& operator=(const CscMatrix<ValueType,IndicesType,base> &C)
        {
            this->setSymmetry(C.isSym());

            this->m_ = C.nRows();
            this->n_ = C.nCols();
            this->nz_ = C.nnz();
            if (this->val_ != 0)
                delete[] this->val_;
            if (this->rows_ != 0)
                delete[] this->rows_;
            if (this->cols_ != 0)
                delete[] this->cols_;

            this->val_ = C.pVal();
            this->rows_ = C.pRows();
            this->cols_ = C.pCols();

            return *this;
        }

        //! Slow get (read-only) of the element at the i-th row and j-th column
        /*!
         * @param i the row indice
         * @param j the column indice
         * @return the element at the i-th row and j-th column
         */
        ValueType get(IndicesType i, IndicesType j) const
        {

            for (IndicesType t = this->cols_[j]; t < this->cols_[j+1]; t++)
                if (this->rows_[t] == i) 
                    return this->val_[t];

            if (this->isSymmetric) {
                for (IndicesType t = this->cols_[i]; t < this->cols_[i+1]; t++)
                    if (this->rows_[t] == j) 
                        return this->val_[t];
            }

            // default behavior if we do not have this entry
            if (i < this->m_ && j < this->n_)
                return 0.0;
            else {
                cerr << "Out of bound access" << endl;
                throw -1;
            }

            return this->val_[0];
        }

        //! Helps to convert without hassle the COO to CSC
        friend void CooToCsc<>(CooMatrix<ValueType,IndicesType,base> &R, CscMatrix<ValueType,IndicesType,base> &L);
};

template<typename ValueType, typename IndicesType, int base>
ostream& operator << (ostream & os, const CscMatrix<ValueType,IndicesType,base> & mat)
{

    IndicesType M = mat.nRows();
    IndicesType N = mat.nCols();
    IndicesType rowp1, colp1;
    int flag = 0;


    std::ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
    std::ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);

    int oldp = os.precision(12);


//  Loop through rows...
    for (IndicesType j = 0; j < N ; j++)
       for (IndicesType i=mat.col(j);j<mat.col(j+1);i++)
       {   
          rowp1 =  mat.row(i) + 1;
          colp1 =  j + 1;
          if ( rowp1 == M && colp1 == N ) flag = 1;
          os.width(14);
          os <<  rowp1 ; os << "    " ;
          os.width(14);
          os <<  colp1 ; os << "    " ;
          os.width(20);
          os <<  mat.val(i) << "\n";
       }

    if (flag == 0)
    {
       os.width(14);
       os <<  M ; os << "    " ;
       os.width(14);
       os <<  N ; os << "    " ;
       os.width(20);
       os <<  mat.get(M-1,N-1) << "\n";
    }

    os.setf(olda,ios::adjustfield);
    os.setf(oldf,ios::floatfield);
    os.precision(oldp);

    return os;
}

#endif //_CSC_HPP_
