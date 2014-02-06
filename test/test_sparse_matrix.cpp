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

TEST (CooMatrix_Constructor, CoordinateConstructor)
{
    int m = 5;
    int n = 5;
    int nz = 8;
    int irn[8] = {0, 1, 2, 2, 3, 3, 4, 4};
    int jcn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
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
    size_t jcn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
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
    int jcn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
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

TEST (CooMatrix_Transpose, CoordinateTranspose)
{
    size_t m = 5;
    size_t n = 5;
    size_t nz = 8;
    size_t irn[8] = {0, 1, 2, 2, 3, 3, 4, 4};
    size_t jcn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CooMatrix<double,size_t,0> cm(m, n, nz, irn, jcn, val);

    CooMatrix<double,size_t,0> cmT = cm.transpose();
    EXPECT_THAT(val, ElementsAreArray(cmT.pVal(),  cm.nnz()));
    EXPECT_THAT(irn, ElementsAreArray(cmT.pCols(), cm.nnz()));
    EXPECT_THAT(jcn, ElementsAreArray(cmT.pRows(), cm.nnz()));

    cm.inPlaceTranspose();
    EXPECT_THAT(val, ElementsAreArray(cm.pVal(),  cm.nnz()));
    EXPECT_THAT(irn, ElementsAreArray(cm.pCols(), cm.nnz()));
    EXPECT_THAT(jcn, ElementsAreArray(cm.pRows(), cm.nnz()));

}

TEST (CooMatrix_Transpose, CoordinateTransposeRef)
{
    size_t m = 5;
    size_t n = 5;
    size_t nz = 8;
    size_t irn[8] = {0, 1, 2, 2, 3, 3, 4, 4};
    size_t jcn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CooMatrix<double,size_t,0> cm(m, n, nz, irn, jcn, val);

    CooMatrix<double,size_t,0> cmTRef = cm.transpose(true);
    EXPECT_THAT(cmTRef.pVal(),Eq(cm.pVal()));
    EXPECT_THAT(cmTRef.pRows(),Eq(cm.pCols()));
    EXPECT_THAT(cmTRef.pCols(),Eq(cm.pRows()));

    cm.inPlaceTranspose();
    EXPECT_THAT(val, ElementsAreArray(cm.pVal(),  cm.nnz()));
    EXPECT_THAT(irn, ElementsAreArray(cm.pCols(), cm.nnz()));
    EXPECT_THAT(jcn, ElementsAreArray(cm.pRows(), cm.nnz()));

}

TEST (CooMatrix_Transpose, CoordinateTransposeInPlace)
{
    size_t m = 5;
    size_t n = 5;
    size_t nz = 8;
    size_t irn[8] = {0, 1, 2, 2, 3, 3, 4, 4};
    size_t jcn[8] = {0, 1, 2, 4, 3, 4, 2, 3};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CooMatrix<double,size_t,0> cm(m, n, nz, irn, jcn, val);

    cm.inPlaceTranspose();
    EXPECT_THAT(val, ElementsAreArray(cm.pVal(),  cm.nnz()));
    EXPECT_THAT(irn, ElementsAreArray(cm.pCols(), cm.nnz()));
    EXPECT_THAT(jcn, ElementsAreArray(cm.pRows(), cm.nnz()));

}

TEST (CooMatrix_Desymmetrize, CoordinateDesymmetrize)
{
    size_t m = 5;
    size_t n = 5;
    size_t nz = 6;
    size_t irn[6] = {0, 1, 2, 3, 4, 4};
    size_t jcn[6] = {0, 1, 2, 3, 2, 3};
    double val[6] = {1, 1, 1, 1, 3, 2};
    
    size_t u_irn[8] = {0, 1, 2, 3, 4, 2, 4, 3};
    size_t u_jcn[8] = {0, 1, 2, 3, 2, 4, 3, 4};
    double u_val[8] = {1, 1, 1, 1, 3, 3, 2, 2};

    CooMatrix<double,size_t,0> cm(m, n, nz, irn, jcn, val);
    cm.setSymmetric();

    EXPECT_THAT(cm.isSym(), Eq(true));
    EXPECT_THAT(cm.get(4,3), Eq(cm.get(3,4)));
    EXPECT_THAT(cm.get(4,4), Eq(0));

    CooMatrix<double,size_t,0> cmDS = cm.desym();
    EXPECT_THAT(cmDS.isSym(), Eq(false));

    EXPECT_THAT(u_val, ElementsAreArray(cmDS.pVal(),  cmDS.nnz()));
    EXPECT_THAT(u_irn, ElementsAreArray(cmDS.pRows(), cmDS.nnz()));
    EXPECT_THAT(u_jcn, ElementsAreArray(cmDS.pCols(), cmDS.nnz()));
}

TEST (CooMatrix_Desymmetrize, CoordinateDesymmetrizeInPlace)
{
    size_t m = 5;
    size_t n = 5;
    size_t nz = 6;
    size_t irn[6] = {0, 1, 2, 3, 4, 4};
    size_t jcn[6] = {0, 1, 2, 3, 2, 3};
    double val[6] = {1, 1, 1, 1, 3, 2};
    
    size_t u_irn[8] = {0, 1, 2, 3, 4, 2, 4, 3};
    size_t u_jcn[8] = {0, 1, 2, 3, 2, 4, 3, 4};
    double u_val[8] = {1, 1, 1, 1, 3, 3, 2, 2};

    CooMatrix<double,size_t,0> cm(m, n, nz, irn, jcn, val);
    cm.setSymmetric();

    EXPECT_THAT(cm.isSym(), Eq(true));
    EXPECT_THAT(cm.get(4,3), Eq(cm.get(3,4)));
    EXPECT_THAT(cm.get(4,4), Eq(0));

    cm.inPlaceDesym();
    EXPECT_THAT(cm.isSym(), Eq(false));

    EXPECT_THAT(u_val, ElementsAreArray(cm.pVal(),  cm.nnz()));
    EXPECT_THAT(u_irn, ElementsAreArray(cm.pRows(), cm.nnz()));
    EXPECT_THAT(u_jcn, ElementsAreArray(cm.pCols(), cm.nnz()));
}

TEST (CooMatrix_Transform, CoordinateToCSR)
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

    CsrMatrix<> rm = cm.toCSR();
    EXPECT_THAT(rm.isSym(),  Eq(true));
    EXPECT_THAT(rm.get(4,3), Eq(rm.get(3,4)));
    EXPECT_THAT(rm.get(4,4), Eq(0));

    EXPECT_THAT(csr_val,  ElementsAreArray(rm.pVal(),  rm.nnz()));
    EXPECT_THAT(csr_iptr, ElementsAreArray(rm.pRows(), rm.nRows() + 1));
    EXPECT_THAT(csr_jcn,  ElementsAreArray(rm.pCols(), rm.nnz()));
}
