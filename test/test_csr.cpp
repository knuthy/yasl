#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "sparse.hpp"
#include "coo.hpp"
#include "csr.hpp"

#include <iostream>

using ::testing::AtLeast;
using ::testing::Return;
using ::testing::Eq;
using ::testing::ElementsAreArray;

TEST (CsrMatrix_Constructor, CSRConstructor)
{
    int m = 5;
    int n = 5;
    int nz = 8;
    int iptr[6] = {0, 1, 2, 4, 6, 8};
    int jcn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CsrMatrix<> cm(m, n, nz, iptr, jcn, val);
    EXPECT_THAT(cm.nRows(), Eq(m));
    EXPECT_THAT(cm.nCols(), Eq(n));
    EXPECT_THAT(cm.nnz(), Eq(nz));
    EXPECT_THAT(val, ElementsAreArray(cm.pVal(),  cm.nnz()));
    EXPECT_THAT(iptr, ElementsAreArray(cm.pRows(), cm.nRows() + 1));
    EXPECT_THAT(jcn, ElementsAreArray(cm.pCols(), cm.nnz()));
    EXPECT_THAT(cm.row(cm.nRows()), Eq(cm.nnz()));
}

TEST (CsrMatrix_Constructor, CSRConstructorCustomType)
{
    size_t m = 5;
    size_t n = 5;
    size_t nz = 8;
    size_t iptr[6] = {0, 1, 2, 4, 6, 8};
    size_t jcn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CsrMatrix<double,size_t,0> cm(m, n, nz, iptr, jcn, val);

    EXPECT_THAT(cm.nRows(), Eq(m));
    EXPECT_THAT(cm.nCols(), Eq(n));
    EXPECT_THAT(cm.nnz(), Eq(nz));
    EXPECT_THAT(val, ElementsAreArray(cm.pVal(),  cm.nnz()));
    EXPECT_THAT(iptr, ElementsAreArray(cm.pRows(), cm.nRows() + 1));
    EXPECT_THAT(jcn, ElementsAreArray(cm.pCols(), cm.nnz()));
    EXPECT_THAT(cm.row(cm.nRows()), Eq(cm.nnz()));
}

TEST (CsrMatrix_Constructor, CSRConstructorReference)
{
    int m = 5;
    int n = 5;
    int nz = 8;
    int iptr[6] = {0, 1, 2, 4, 6, 8};
    int jcn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CsrMatrix<> cm(m, n, nz, iptr, jcn, val, true);

    //samething as before, except be sure that
    //the arrays are pointing to the same memory location
    EXPECT_THAT(cm.nRows(), Eq(m));
    EXPECT_THAT(cm.nCols(), Eq(n));
    EXPECT_THAT(cm.nnz(), Eq(nz));
    EXPECT_THAT(val, ElementsAreArray(cm.pVal(),  cm.nnz()));
    EXPECT_THAT(iptr, ElementsAreArray(cm.pRows(), cm.nRows() + 1));
    EXPECT_THAT(jcn, ElementsAreArray(cm.pCols(), cm.nnz()));
    EXPECT_THAT(cm.pVal(),Eq(val));
    EXPECT_THAT(cm.pRows(),Eq(iptr));
    EXPECT_THAT(cm.pCols(),Eq(jcn));
}

TEST (CsrMatrix_Transform, CSRfromCoordinate)
{
    int m = 5;
    int n = 5;
    int nz = 8;

    int irn[8] = {0, 1, 2, 2, 3, 3, 4, 4};
    int jcn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    int csr_iptr[6] = {0, 1, 2, 4, 6, 8};
    int csr_jcn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
    double csr_val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CooMatrix<> cm(m, n, nz, irn, jcn, val);
    cm.setSymmetric();
    EXPECT_THAT(cm.isSym(), Eq(true));
    EXPECT_THAT(cm.get(4,3), Eq(cm.get(3,4)));
    EXPECT_THAT(cm.get(4,4), Eq(0));

    CsrMatrix<> vm(cm);

    CsrMatrix<> rm = cm.toCSR();
    EXPECT_THAT(rm.isSym(),  Eq(true));
    EXPECT_THAT(rm.get(4,3), Eq(rm.get(3,4)));
    EXPECT_THAT(rm.get(4,4), Eq(0));

    EXPECT_THAT(csr_val,  ElementsAreArray(rm.pVal(),  rm.nnz()));
    EXPECT_THAT(csr_iptr, ElementsAreArray(rm.pRows(), rm.nRows() + 1));
    EXPECT_THAT(csr_jcn,  ElementsAreArray(rm.pCols(), rm.nnz()));
}
