#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "sparse.hpp"
#include "coo.hpp"
#include "csr.hpp"
#include "csc.hpp"

#include <iostream>

using ::testing::AtLeast;
using ::testing::Return;
using ::testing::Eq;
using ::testing::ElementsAreArray;

TEST (CscMatrix_Constructor, CSCConstructor)
{
    int m = 5;
    int n = 5;
    int nz = 8;
    int irn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
    int jptr[6] = {0, 1, 2, 4, 6, 8};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CscMatrix<> cm(m, n, nz, irn, jptr, val);
    EXPECT_THAT(cm.nRows(), Eq(m));
    EXPECT_THAT(cm.nCols(), Eq(n));
    EXPECT_THAT(cm.nnz(), Eq(nz));
    EXPECT_THAT(val, ElementsAreArray(cm.pVal(),  cm.nnz()));
    EXPECT_THAT(irn, ElementsAreArray(cm.pRows(), cm.nnz()));
    EXPECT_THAT(jptr, ElementsAreArray(cm.pCols(), cm.nCols() + 1));
    EXPECT_THAT(cm.col(cm.nCols()), Eq(cm.nnz()));
}

TEST (CscMatrix_Constructor, CSCConstructorCustomType)
{
    size_t m = 5;
    size_t n = 5;
    size_t nz = 8;
    size_t jptr[6] = {0, 1, 2, 4, 6, 8};
    size_t irn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CscMatrix<double,size_t,0> cm(m, n, nz, irn, jptr, val);

    EXPECT_THAT(cm.nRows(), Eq(m));
    EXPECT_THAT(cm.nCols(), Eq(n));
    EXPECT_THAT(cm.nnz(), Eq(nz));
    EXPECT_THAT(val, ElementsAreArray(cm.pVal(),  cm.nnz()));
    EXPECT_THAT(irn, ElementsAreArray(cm.pRows(), cm.nnz()));
    EXPECT_THAT(jptr, ElementsAreArray(cm.pCols(), cm.nCols() + 1));
    EXPECT_THAT(cm.col(cm.nCols()), Eq(cm.nnz()));
}

TEST (CscMatrix_Constructor, CSCConstructorReference)
{
    int m = 5;
    int n = 5;
    int nz = 8;
    int jptr[6] = {0, 1, 2, 4, 6, 8};
    int irn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CscMatrix<> cm(m, n, nz, irn, jptr, val, true);

    //samething as before, except be sure that
    //the arrays are pointing to the same memory location
    EXPECT_THAT(cm.nRows(), Eq(m));
    EXPECT_THAT(cm.nCols(), Eq(n));
    EXPECT_THAT(cm.nnz(), Eq(nz));
    EXPECT_THAT(val, ElementsAreArray(cm.pVal(),  cm.nnz()));
    EXPECT_THAT(irn, ElementsAreArray(cm.pRows(), cm.nnz()));
    EXPECT_THAT(jptr, ElementsAreArray(cm.pCols(), cm.nCols() + 1));
    EXPECT_THAT(cm.pVal(),Eq(val));
    EXPECT_THAT(cm.pRows(),Eq(irn));
    EXPECT_THAT(cm.pCols(),Eq(jptr));
}

TEST (CscMatrix_Transform, CSCfromCoordinate)
{
    int m = 5;
    int n = 5;
    int nz = 8;

    int irn[8] = {0, 1, 2, 2, 3, 3, 4, 4};
    int jcn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CooMatrix<> cm(m, n, nz, irn, jcn, val);
    cm.setSymmetric();
    EXPECT_THAT(cm.isSym(), Eq(true));
    EXPECT_THAT(cm.get(4,3), Eq(cm.get(3,4)));
    EXPECT_THAT(cm.get(4,4), Eq(0));

    CscMatrix<> rm = cm.toCSC();

    EXPECT_THAT(rm.isSym(),  Eq(true));
    EXPECT_THAT(rm.get(4,3), Eq(rm.get(3,4)));
    EXPECT_THAT(rm.get(4,4), Eq(0));

    int csc_irn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
    int csc_jptr[6] = {0, 1, 2, 4, 6, 8};
    double csc_val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    EXPECT_THAT(csc_val,  ElementsAreArray(rm.pVal(),  rm.nnz()));
    EXPECT_THAT(csc_irn,  ElementsAreArray(rm.pRows(), rm.nnz()));
    EXPECT_THAT(csc_jptr, ElementsAreArray(rm.pCols(), rm.nCols() + 1));
}
