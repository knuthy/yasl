#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "sparse.hpp"
#include "coo.hpp"

#include <iostream>

using ::testing::AtLeast;
using ::testing::Return;
using ::testing::Eq;
using ::testing::ElementsAreArray;

TEST (CooMatrix_Constructor, CoordinateConstructor)
{
    int m = 5;
    int n = 5;
    int nz = 8;
    int irn[8] = {0, 1, 2, 2, 3, 3, 4, 4};
    int jcn[8] = {0, 1, 2, 4, 3, 4, 2, 1};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CooMatrix<> cm(m, n, nz, irn, jcn, val);
    EXPECT_THAT(cm.nRows(), Eq(m));
    EXPECT_THAT(cm.nCols(), Eq(n));
    EXPECT_THAT(cm.nnz(), Eq(nz));
    EXPECT_THAT(val, ElementsAreArray(cm.pVal(),  cm.nnz()));
    EXPECT_THAT(irn, ElementsAreArray(cm.pRows(), cm.nnz()));
    EXPECT_THAT(jcn, ElementsAreArray(cm.pCols(), cm.nnz()));
}

TEST (CooMatrix_Constructor, CoordinateConstructorCustomType)
{
    size_t m = 5;
    size_t n = 5;
    size_t nz = 8;
    size_t irn[8] = {0, 1, 2, 2, 3, 3, 4, 4};
    size_t jcn[8] = {0, 1, 2, 4, 3, 4, 2, 1};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CooMatrix<double,size_t,0> cm(m, n, nz, irn, jcn, val);
    EXPECT_THAT(cm.nRows(), Eq(m));
    EXPECT_THAT(cm.nCols(), Eq(n));
    EXPECT_THAT(cm.nnz(), Eq(nz));
    EXPECT_THAT(val, ElementsAreArray(cm.pVal(),  cm.nnz()));
    EXPECT_THAT(irn, ElementsAreArray(cm.pRows(), cm.nnz()));
    EXPECT_THAT(jcn, ElementsAreArray(cm.pCols(), cm.nnz()));
}

TEST (CooMatrix_Constructor, CoordinateConstructorReference)
{
    int m = 5;
    int n = 5;
    int nz = 8;
    int irn[8] = {0, 1, 2, 2, 3, 3, 4, 4};
    int jcn[8] = {0, 1, 2, 4, 3, 4, 2, 1};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CooMatrix<> cm(m, n, nz, irn, jcn, val, true);

    //samething as before, except be sure that
    //the arrays are pointing to the same memory location
    EXPECT_THAT(cm.nRows(), Eq(m));
    EXPECT_THAT(cm.nCols(), Eq(n));
    EXPECT_THAT(cm.nnz(), Eq(nz));
    EXPECT_THAT(val, ElementsAreArray(cm.pVal(),  cm.nnz()));
    EXPECT_THAT(irn, ElementsAreArray(cm.pRows(), cm.nnz()));
    EXPECT_THAT(jcn, ElementsAreArray(cm.pCols(), cm.nnz()));
    EXPECT_THAT(cm.pVal(),Eq(val));
    EXPECT_THAT(cm.pRows(),Eq(irn));
    EXPECT_THAT(cm.pCols(),Eq(jcn));
}

TEST (CooMatrix_Constructor, CoordinateTranspose)
{
    size_t m = 5;
    size_t n = 5;
    size_t nz = 8;
    size_t irn[8] = {0, 1, 2, 2, 3, 3, 4, 4};
    size_t jcn[8] = {0, 1, 2, 4, 3, 4, 2, 1};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CooMatrix<double,size_t,0> cm(m, n, nz, irn, jcn, val);

    CooMatrix<double,size_t,0> cmTRef = cm.transpose(true);
    EXPECT_THAT(cmTRef.pVal(),Eq(cm.pVal()));
    EXPECT_THAT(cmTRef.pRows(),Eq(cm.pCols()));
    EXPECT_THAT(cmTRef.pCols(),Eq(cm.pRows()));

    CooMatrix<double,size_t,0> cmT = cm.transpose();
    EXPECT_THAT(val, ElementsAreArray(cmT.pVal(),  cm.nnz()));
    EXPECT_THAT(irn, ElementsAreArray(cmT.pCols(), cm.nnz()));
    EXPECT_THAT(jcn, ElementsAreArray(cmT.pRows(), cm.nnz()));
}

