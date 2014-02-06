#ifndef _COO_HPP_
#define _COO_HPP_

#include <vector>
#include <iostream>
#include <stdlib.h>

#include "sparse.hpp"

using namespace std;

template< typename ValueType, typename IndicesType, int base>
void CooToCsr(CooMatrix<ValueType,IndicesType,base> &R, CsrMatrix<ValueType,IndicesType,base> &L);

template< typename ValueType, typename IndicesType, int base>
void CooToCsc(CooMatrix<ValueType,IndicesType,base> &R, CscMatrix<ValueType,IndicesType,base> &L);

//! Coordinate storage format sparse matrix
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
class CooMatrix : public SparseMatrix<ValueType,IndicesType,base> {
    protected:
        //! Creates an coordinate storage sparse matrix with vectors as reference
        /*!
         * @param m number of rows
         * @param n number of columns
         * @param nz number of nonzero entries
         * @param irn the row indices vector
         * @param jcn the column indices vector
         * @param val the entries vector
         * @param reference the vectors in the matrix are the real ones 
         * @param switchToBase if @base == 0 then apply a -1 to indices, if @base == 1 then apply a +1
         */
        void init( IndicesType m,  IndicesType n,  IndicesType nz,
                 IndicesType* irn,  IndicesType* jcn,  ValueType* val,
                 bool reference, bool swithToBase)
        {
            this->m_ = m;
            this->n_ = n;
            this->nz_ = nz;
            this->isReference = reference;
            if (this->isReference) {
                this->val_  = val;
                this->rows_ = irn;
                this->cols_ = jcn;
            } else {
                this->val_ = new ValueType[this->nz_];
                this->rows_ = new IndicesType[this->nz_];
                this->cols_ = new IndicesType[this->nz_];
                if (swithToBase) {
                    if (base == 0)
                        for (IndicesType t = 0; t < this->nz_; t++) {
                            this->val_[t] = val[t];
                            this->rows_[t] = irn[t] - 1;
                            this->cols_[t] = jcn[t] - 1;
                        }
                    else
                        for (IndicesType t = 0; t < this->nz_; t++) {
                            this->val_[t] = val[t];
                            this->rows_[t] = irn[t] + 1;
                            this->cols_[t] = jcn[t] + 1;
                        }
                } else {
                    std::copy(val, val + this->nz_, this->val_);
                    std::copy(irn, irn + this->nz_, this->rows_);
                    std::copy(jcn, jcn + this->nz_, this->cols_);
                }
            }
            this->isSymmetric = false;
        }
    public:
        /*!
         * Basic constructor, initialize everything to zero
         */
        CooMatrix()
        {
            this->m_ = 0;
            this->n_ = 0;
            this->nz_ = 0;
            this->val_ = NULL;
            this->rows_ = NULL;
            this->cols_ = NULL;
            this->isReference = false;
            this->isSymmetric = false;
        }

        //! Creates an empty coordinate storage sparse matrix
        /*!
         * @param m number of rows
         * @param n number of columns
         * @param nz number of nonzero entries
         */
        CooMatrix( IndicesType m,  IndicesType n,  IndicesType nz)
        {
            this->m_ = m;
            this->n_ = n;
            this->nz_ = nz;
            this->val_ = new ValueType[this->nz_];
            this->rows_ = new IndicesType[this->nz_];
            this->cols_ = new IndicesType[this->nz_];
            this->isReference = false;
            this->isSymmetric = false;
        }


        //! Creates an coordinate storage sparse matrix
        /*!
         * @param m number of rows
         * @param n number of columns
         * @param nz number of nonzero entries
         * @param irn the row indices vector
         * @param jcn the column indices vector
         * @param val the entries vector
         */
        CooMatrix( IndicesType m,  IndicesType n,  IndicesType nz,
                 IndicesType* irn,  IndicesType* jcn,  ValueType* val)
        {
            this->init(m, n, nz, irn, jcn, val, false, false);
        }

        //! Creates an coordinate storage sparse matrix
        /*!
         * @param m number of rows
         * @param n number of columns
         * @param nz number of nonzero entries
         * @param irn the row indices vector
         * @param jcn the column indices vector
         * @param val the entries vector
         * @param reference the vectors in the matrix are the real ones 
         */
        CooMatrix( IndicesType m,  IndicesType n,  IndicesType nz,
                 IndicesType* irn,  IndicesType* jcn,  ValueType* val,
                 bool reference)
        {
            this->init(m, n, nz, irn, jcn, val, reference, false);
        }

        //! Creates an coordinate storage sparse matrix
        /*!
         * @param m number of rows
         * @param n number of columns
         * @param nz number of nonzero entries
         * @param irn the row indices vector
         * @param jcn the column indices vector
         * @param val the entries vector
         * @param reference the vectors in the matrix are the real ones 
         * @param switchToBase if @base == 0 then apply a -1 to indices, if @base == 1 then apply a +1
         */
        CooMatrix( IndicesType m,  IndicesType n,  IndicesType nz,
                 IndicesType* irn,  IndicesType* jcn,  ValueType* val,
                 bool reference, bool swithToBase)
        {
            this->init( m,  n,  nz, irn, jcn, val, reference, swithToBase);
        }

        //! Assign the right-hand side matrix into the left-hand side (reference)
        /*!
         * Returns a reference to the left-hand side matrix, if you want a copy
         * use the constructor. If the matrix contained data, it will be deleted.
         *
         * @return a reference to the left-hand side
         */
        CooMatrix& operator=(const CooMatrix &C)
        {
            this->setSymmetry(C.isSym());

            this->m_ = C.nRows();
            this->n_ = C.nCols();

            if (this->nz_ != C.nnz()) {
                this->nz_ = C.nnz();
                if (this->val_ != NULL)
                    delete[] this->val_;
                if (this->rows_ != NULL)
                    delete[] this->rows_;
                if (this->cols_ != NULL)
                    delete[] this->cols_;

                this->val_ = new ValueType[this->nz_];
                this->rows_ = new IndicesType[this->nz_];
                this->cols_ = new IndicesType[this->nz_];
            } else {
                // be sure that everything is allocated
                if (this->val_ == NULL)
                    this->val_ = new ValueType[this->nz_];
                if (this->rows_ == NULL)
                    this->rows_ = new IndicesType[this->nz_];
                if (this->cols_ == NULL)
                    this->cols_ = new IndicesType[this->nz_];
            }

            std::copy(C.pVal(), C.pVal() + this->nz_, this->val_);
            std::copy(C.pRows(), C.pRows() + this->nz_, this->rows_);
            std::copy(C.pCols(), C.pCols() + this->nz_, this->cols_);

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

            if (this->isSymmetric) {
                for (IndicesType t = 0; t < this->nz_; t++)
                    if ((this->rows_[t] == i && this->cols_[t] == j) ||
                        (this->rows_[t] == j && this->cols_[t] == i))
                        return this->val_[t];

            } else {
                for (IndicesType t = 0; t < this->nz_; t++)
                    if (this->rows_[t] == i && this->cols_[t] == j) 
                        return this->val_[t];
            }

            // default behavior if we do not have this entry
            if (i < this->m_ && j < this->n_) return 0.0;

            return this->val_[0];
        }

        //! Slow access (read/write) to the element at the i-th row and j-th column
        /*!
         * @param i the row indice
         * @param j the column indice
         * @return the element at the i-th row and j-th column
         */
        ValueType& operator() (IndicesType i, IndicesType j)
        {
            for (IndicesType t = 0; t < this->nz_; t++)
                if (this->rows_[t] == i && this->cols_[t] == j) return this->val_[t];

            std::cerr << "Array element (" << i << "," << j ;
            std::cerr << ") not in sparse structure -- cannot access." << "\n";
            /// @todo: define exceptions
            exit(1);
            return this->val_[0];
        }

        //! Transpose the matrix
        /*!
         * @param refrence if true, the vectors in the resulting matrix share the same memory
         * @return the transpose of the matrix
         */
        CooMatrix transpose(bool reference) const
        {
            return CooMatrix<ValueType,IndicesType,base>(this->n_, this->m_, this->nz_,
                    this->cols_, this->rows_, this->val_, reference);
        }

        //! Transpose the matrix
        /*!
         * @return a copy of the transpose of the matrix
         */
        CooMatrix transpose() const
        {
            return this->transpose(false);
        }

        //! Transpose the matrix, inplace
        /*!
         * @return the transpose of the matrix
         */
        CooMatrix& inPlaceTranspose(){
            *this = this->transpose();
            return *this;
        }

        //! Desymmetrize the matrix
        /*!
         * @return the matrix in unsymmetric format
         */
        CooMatrix desym(){
            if (! this->isSymmetric) {
                cerr << "Error: the matrix is already non-symmetric" << endl;
                throw -1; // @todo should have a better exception handling
            }

            IndicesType *t_irn = new IndicesType[2*this->nz_];
            IndicesType *t_jcn = new IndicesType[2*this->nz_];
            ValueType   *t_val = new ValueType[2*this->nz_];
            IndicesType t_nz = 0;

            for(IndicesType k = 0; k < this->nz_; ++k) {
                t_irn[t_nz] = this->rows_[k];
                t_jcn[t_nz] = this->cols_[k];
                t_val[t_nz] = this->val_[k];
                t_nz++;
                if(this->rows_[k] != this->cols_[k]){
                    t_irn[t_nz] = this->cols_[k];
                    t_jcn[t_nz] = this->rows_[k];
                    t_val[t_nz] = this->val_[k];
                    t_nz++;
                }
            }
            CooMatrix<ValueType,IndicesType,base> ds(
                    this->m_, this->n_, t_nz,
                    t_irn, t_jcn, t_val);
            delete[] t_irn;
            delete[] t_jcn;
            delete[] t_val;
            return ds;
        }

        //! Desymmetrize the matrix, in place
        /*!
         * @return the matrix in unsymmetric format
         */
        CooMatrix& inPlaceDesym(){
            *this = this->desym();
            return *this;
        }

        //! Converts the matrix to compressed row storage format
        /*!
         * @return the matrix in compressed row stroage format
         */
        CsrMatrix<ValueType,IndicesType,base> toCSR()
        {
            CsrMatrix<ValueType,IndicesType,base> csr(this->m_, this->n_, this->nz_);
            CooToCsr(*this, csr);
            return csr;
        }

        //! Converts the matrix to compressed column storage format
        /*!
         * @return the matrix in compressed column stroage format
         */
        CscMatrix<ValueType,IndicesType,base> toCSC()
        {
            CscMatrix<ValueType,IndicesType,base> csc(this->m_, this->n_, this->nz_);
            CooToCsc(*this, csc);
            return csc;
        }

        //! Helps to convert without hassle the COO to CSR
        friend void CooToCsr<>(CooMatrix<ValueType,IndicesType,base> &R, CsrMatrix<ValueType,IndicesType,base> &L);
        friend void CooToCsc<>(CooMatrix<ValueType,IndicesType,base> &R, CscMatrix<ValueType,IndicesType,base> &L);
};

template<typename ValueType, typename IndicesType, int base>
ostream& operator << (ostream & os, const CooMatrix<ValueType,IndicesType,base> & mat)
{
         IndicesType nnz = mat.nnz();
         IndicesType M = mat.nRows();
         IndicesType N = mat.nCols();
         IndicesType rowp1, colp1;
         int flag = 0;


        std::ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
        std::ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);

         int oldp = os.precision(12);

//       Loop through Nonzeros
         for (IndicesType j = 0; j < nnz ; j++) 
         {
            rowp1 = mat.row(j) +1;
            colp1 = mat.col(j) +1;
            if ( rowp1 == M && colp1 == N ) flag = 1;
            os.width(14);
            os <<  rowp1 ; os << "    " ;
            os.width(14);
            os <<  colp1 ; os << "    " ;
            os.width(20);
            os <<  mat.val(j) << "\n";
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

template< typename ValueType, typename IndicesType, int base>
void CooToCsr(CooMatrix<ValueType,IndicesType,base> &R, CsrMatrix<ValueType,IndicesType,base> &L)
{
    // same behaviour as scipy coo->csr 
    // (see scipy/sparse/sparsetools/coo.h)
    IndicesType *Bp = L.pRows();
    IndicesType *Bj = L.pCols();
    ValueType *Bx   = L.pVal();

    //compute number of non-zero entries per row of A
    std::fill(Bp, Bp + R.m_, 0);

    for (IndicesType t = 0; t < R.nz_; t++){
        Bp[R.rows_[t]]++;
    }

    //cumsum the nnz per row to get Bp[]
    for(IndicesType i = 0, cumsum = 0; i < R.m_; i++){
        IndicesType temp = Bp[i];
        Bp[i] = cumsum;
        cumsum += temp;
    }
    Bp[R.m_] = R.nz_;

    for(IndicesType t = 0; t < R.nz_; t++){
        IndicesType row = R.rows_[t];
        IndicesType dest = Bp[row];

        Bj[dest] = R.cols_[t];
        Bx[dest] = R.val_[t];

        Bp[row]++;
    }

    for(IndicesType i = 0, last = 0; i <= R.m_; i++){
        IndicesType temp = Bp[i];
        Bp[i] = last;
        last = temp;
    }

    L.setSymmetry(R.isSymmetric);
}
template< typename ValueType, typename IndicesType, int base>
void CooToCsc(CooMatrix<ValueType,IndicesType,base> &R, CscMatrix<ValueType,IndicesType,base> &L)
{
    // same behaviour as scipy coo->csc 
    // (see scipy/sparse/sparsetools/coo.h)
    IndicesType *Bp = L.pCols();
    IndicesType *Bi = L.pRows();
    ValueType *Bx   = L.pVal();

    //compute number of non-zero entries per row of A
    std::fill(Bp, Bp + R.n_, 0);

    for (IndicesType t = 0; t < R.nz_; t++){
        Bp[R.cols_[t]]++;
    }

    //cumsum the nnz per row to get Bp[]
    for(IndicesType i = 0, cumsum = 0; i < R.n_; i++){
        IndicesType temp = Bp[i];
        Bp[i] = cumsum;
        cumsum += temp;
    }
    Bp[R.n_] = R.nz_;

    for(IndicesType t = 0; t < R.nz_; t++){
        IndicesType col = R.cols_[t];
        IndicesType dest = Bp[col];

        Bi[dest] = R.rows_[t];
        Bx[dest] = R.val_[t];

        Bp[col]++;
    }

    for(IndicesType i = 0, last = 0; i <= R.n_; i++){
        IndicesType temp = Bp[i];
        Bp[i] = last;
        last = temp;
    }

    L.setSymmetry(R.isSymmetric);
}

#endif // _COO_HPP_
