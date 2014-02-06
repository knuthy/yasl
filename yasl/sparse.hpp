#ifndef _SPARSE_HPP_
#define _SPARSE_HPP_

#include <vector>
#include <iostream>

//! Sparse matrix class
/*!
 * @tparam ValueType the type of the entries of the matrix, for the moment only double is supported
 * @tparam IndicesType the type of vector that stores row/column indices, default is size_t
 * @tparam base the base of the vector indices, 0 for c-based, 1 for fortran-based
 */
template< typename ValueType, typename IndicesType, int base = 0>
class SparseMatrix {
    protected:
        ValueType* val_; /*!< data values (nz_ elements)*/
        IndicesType* rows_; /*!< row indices or pointers, depends on the storage format (m_ or nz_ elements) */
        IndicesType* cols_; /*!< column indices or pointers, depends on the storage format (n_ or nz_ elements) */

        IndicesType nz_; /*!< the number of nonzero elements */
        IndicesType m_; /*!< the number of rows */
        IndicesType n_; /*!< the number of columns */

        bool isSymmetric; /*!< defines wether the matrix is symmetric or not */
        bool isReference; /*!< defines if the matrix is a reference to dynamic arrays */

    public:
        /*!
         * Basic constructor, initialize everything to zero
         */
        SparseMatrix() : isSymmetric(false), isReference(false), m_(0), n_(0), nz_(0) { }

        ~SparseMatrix()
        {
            // avoid double free problems by not freeing objects that belongs to others
            if (!isReference) {
                delete[] val_;
                delete[] rows_;
                delete[] cols_;
            }
        }

        //! Returns a pointer to the first element of the values vector
        ValueType *pVal() const { return val_; }
        //! Returns a pointer to the first element of the rows vector
        IndicesType *pRows() const { return rows_; }
        //! Returns a pointer to the first element of the cols vector
        IndicesType *pCols() const { return cols_; }

        //! Return the number of rows
        IndicesType nRows() const { return this->m_; };
        //! Return the number of columns
        IndicesType nCols() const { return this->n_; };
        //! Return the number of nonzero elements
        IndicesType nnz() const { return this->nz_; };

        //! Define the symmetry of the matrix
        void setSymmetry(bool isSym) { this->isSymmetric = isSym; }
        //! Define the matrix as symmetric
        void setSymmetric() { this->isSymmetric = true; }
        //! Define the matrix as non-symmetric
        void setNonSymmetric() { this->isSymmetric = false; }
        //! Get the symmetry
        bool isSym() const { return this->isSymmetric; }

        //! Get (read-only) the element at the i-th row and j-th column
        /*!
         * @param i the row indice
         * @param j the column indice
         * @return the element at the i-th row and j-th column
         */
        ValueType get(IndicesType i, IndicesType j) const;

        //! Access (read/write) the element at the i-th row and j-th column
        /*!
         * @param i the row indice
         * @param j the column indice
         * @return the element at the i-th row and j-th column
         */
        ValueType& operator() (IndicesType i, IndicesType j) const;


        //! Returns the entry at position @i
        ValueType val(IndicesType i) const { return this->val_[i]; }
        //! Returns the row indice at position @i
        IndicesType row(IndicesType i) const { return this->rows_[i]; }
        //! Returns the column indice at position @i
        IndicesType col(IndicesType i) const { return this->cols_[i]; }
};

template<typename ValueType, typename IndicesType,int base >
class CooMatrix;

template<typename ValueType, typename IndicesType,int base >
class CscMatrix;

template<typename ValueType, typename IndicesType,int base >
class CsrMatrix;


#endif // _SPARSE_HPP_
