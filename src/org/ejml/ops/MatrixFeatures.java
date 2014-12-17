/*
 * Copyright (c) 2009-2014, Peter Abeles. All Rights Reserved.
 *
 * This file is part of Efficient Java Matrix Library (EJML).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.ejml.ops;

import org.ejml.UtilEjml;
import org.ejml.alg.dense.decomposition.chol.CholeskyDecompositionInner_D32;
import org.ejml.alg.dense.mult.VectorVectorMult;
import org.ejml.data.Complex32F;
import org.ejml.data.D1Matrix32F;
import org.ejml.data.DenseMatrix32F;
import org.ejml.data.ReshapeMatrix32F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.interfaces.decomposition.LUDecomposition;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;


/**
 * <p>
 * Used to compute features that describe the structure of a matrix.
 * <p>
 *
 * <p>
 * Unless explicitly stated otherwise it is assumed that the elements of input matrices
 * contain only real numbers.  If an element is NaN or infinite then the behavior is undefined.
 * See IEEE 754 for more information on this issue.
 * </p>
 *
 * @author Peter Abeles
 */
public class MatrixFeatures {

    /**
     * Checks to see if any element in the matrix is NaN.
     *
     * @param m A matrix. Not modified.
     * @return True if any element in the matrix is NaN.
     */
    public static boolean hasNaN( D1Matrix32F m )
    {
        int length = m.getNumElements();

        for( int i = 0; i < length; i++ ) {
            if( Float.isNaN(m.get(i)))
                return true;
        }
        return false;
    }

    /**
     * Checks to see if any element in the matrix is NaN of Infinite.
     *
     * @param m A matrix. Not modified.
     * @return True if any element in the matrix is NaN of Infinite.
     */
    public static boolean hasUncountable( D1Matrix32F m )
    {
        int length = m.getNumElements();

        for( int i = 0; i < length; i++ ) {
            float a = m.get(i);
            if( Float.isNaN(a) || Float.isInfinite(a))
                return true;
        }
        return false;
    }

    /**
     * Checks to see all the elements in the matrix are zeros
     *
     * @param m A matrix. Not modified.
     * @return True if all elements are zeros or false if not
     */
    public static boolean isZeros( D1Matrix32F m , float tol )
    {
        int length = m.getNumElements();

        for( int i = 0; i < length; i++ ) {
            if( Math.abs(m.get(i)) > tol )
                return false;
        }
        return true;
    }

    /**
     * Checks to see if the matrix is a vector or not.
     *
     * @param mat A matrix. Not modified.
     *
     * @return True if it is a vector and false if it is not.
     */
    public static boolean isVector( D1Matrix32F mat ) {
        return (mat.numCols == 1 || mat.numRows == 1);
    }

    /**
     * <p>
     * Checks to see if the matrix is positive definite.
     * </p>
     * <p>
     * x<sup>T</sup> A x > 0<br>
     * for all x where x is a non-zero vector and A is a symmetric matrix.
     * </p>
     *
     * @param A square symmetric matrix. Not modified.
     *
     * @return True if it is positive definite and false if it is not.
     */
    public static boolean isPositiveDefinite( DenseMatrix32F A ) {
        if( !isSquare(A))
           return false;

        CholeskyDecompositionInner_D32 chol = new CholeskyDecompositionInner_D32(true);
        if( chol.inputModified() )
            A = A.copy();

        return chol.decompose(A);
    }

    /**
     * <p>
     * Checks to see if the matrix is positive semidefinite:
     * </p>
     * <p>
     * x<sup>T</sup> A x >= 0<br>
     * for all x where x is a non-zero vector and A is a symmetric matrix.
     * </p>
     *
     * @param A square symmetric matrix. Not modified.
     *
     * @return True if it is positive semidefinite and false if it is not.
     */
    public static boolean isPositiveSemidefinite( DenseMatrix32F A ) {
        if( !isSquare(A))
           return false;

        EigenDecomposition<DenseMatrix32F> eig = DecompositionFactory.eig(A.numCols,false);
        if( eig.inputModified() )
            A = A.copy();
        eig.decompose(A);

        for( int i = 0; i < A.numRows; i++ ) {
            Complex32F v = eig.getEigenvalue(i);

            if( v.getReal() < 0 )
                return false;
        }

        return true;
    }

    /**
     * Checks to see if it is a square matrix.  A square matrix has
     * the same number of rows and columns.
     *
     * @param mat A matrix. Not modified.
     * @return True if it is a square matrix and false if it is not.
     */
    public static boolean isSquare( D1Matrix32F mat ) {
        return mat.numCols == mat.numRows;
    }

    /**
     * <p>
     * Returns true if the matrix is symmetric within the tolerance.  Only square matrices can be
     * symmetric.
     * </p>
     * <p>
     * A matrix is symmetric if:<br>
     * |a<sub>ij</sub> - a<sub>ji</sub>| &le; tol
     * </p>
     *
     * @param m A matrix. Not modified.
     * @param tol Tolerance for how similar two elements need to be.
     * @return true if it is symmetric and false if it is not.
     */
    public static boolean isSymmetric( DenseMatrix32F m , float tol ) {
        if( m.numCols != m.numRows )
            return false;

        float max = CommonOps.elementMaxAbs(m);

        for( int i = 0; i < m.numRows; i++ ) {
            for( int j = 0; j < i; j++ ) {
                float a = m.get(i,j)/max;
                float b = m.get(j,i)/max;

                float diff = Math.abs(a-b);

                if( !(diff <= tol) ) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * <p>
     * Returns true if the matrix is perfectly symmetric.  Only square matrices can be symmetric.
     * </p>
     * <p>
     * A matrix is symmetric if:<br>
     * a<sub>ij</sub> == a<sub>ji</sub>
     * </p>
     *
     * @param m A matrix. Not modified.
     * @return true if it is symmetric and false if it is not.
     */
    public static boolean isSymmetric( DenseMatrix32F m ) {
        return isSymmetric(m,0.0f);
    }

    /**
     * <p>
     * Checks to see if a matrix is skew symmetric with in tolerance:<br>
     * <br>
     * -A = A<sup>T</sup><br>
     * or<br>
     * |a<sub>ij</sub> + a<sub>ji</sub>| &le; tol
     * </p>
     *
     * @param A The matrix being tested.
     * @param tol Tolerance for being skew symmetric.
     * @return True if it is skew symmetric and false if it is not.
     */
    public static boolean isSkewSymmetric( DenseMatrix32F A , float tol ){
        if( A.numCols != A.numRows )
            return false;

        for( int i = 0; i < A.numRows; i++ ) {
            for( int j = 0; j < i; j++ ) {
                float a = A.get(i,j);
                float b = A.get(j,i);

                float diff = Math.abs(a+b);

                if( !(diff <= tol) ) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Checks to see if the two matrices are inverses of each other.
     *
     * @param a A matrix. Not modified.
     * @param b A matrix. Not modified.
     */
    public static boolean isInverse( DenseMatrix32F a , DenseMatrix32F b , float tol ) {
        if( a.numRows != b.numRows || a.numCols != b.numCols ) {
            return false;
        }

        int numRows = a.numRows;
        int numCols = a.numCols;

        for( int i = 0; i < numRows; i++ ) {
            for( int j = 0; j < numCols; j++ ) {
                float total = 0;
                for( int k = 0; k < numCols; k++ ) {
                    total += a.get(i,k)*b.get(k,j);
                }

                if( i == j ) {
                    if( !(Math.abs(total-1) <= tol) )
                        return false;
                } else if( !(Math.abs(total) <= tol) )
                    return false;
            }
        }

        return true;
    }

    /**
     * <p>
     * Checks to see if each element in the two matrices are within tolerance of
     * each other: tol &ge; |a<sub>ij</sub> - b<sub>ij</sub>|.
     * <p>
     *
     * <p>
     * NOTE: If any of the elements are not countable then false is returned.<br>
     * NOTE: If a tolerance of zero is passed in this is equivalent to calling
     * {@link #isEquals(org.ejml.data.D1Matrix32F, org.ejml.data.D1Matrix32F)}
     * </p>
     *
     * @param a A matrix. Not modified.
     * @param b A matrix. Not modified.
     * @param tol How close to being identical each element needs to be.
     * @return true if equals and false otherwise.
     */
    public static boolean isEquals( D1Matrix32F a , D1Matrix32F b , float tol )
    {
        if( a.numRows != b.numRows || a.numCols != b.numCols ) {
            return false;
        }

        if( tol == 0.0 )
            return isEquals(a,b);

        final int length = a.getNumElements();

        for( int i = 0; i < length; i++ ) {
            if( !(tol >= Math.abs(a.get(i) - b.get(i))) ) {
                return false;
            }
        }
        return true;
    }

    /**
     * <p>
     * Checks to see if each element in the upper or lower triangular portion of the two matrices are within tolerance of
     * each other: tol &ge; |a<sub>ij</sub> - b<sub>ij</sub>|.
     * <p>
     *
     * <p>
     * NOTE: If any of the elements are not countable then false is returned.<br>
     * NOTE: If a tolerance of zero is passed in this is equivalent to calling
     * {@link #isEquals(org.ejml.data.D1Matrix32F, org.ejml.data.D1Matrix32F)}
     * </p>
     *
     * @param a A matrix. Not modified.
     * @param b A matrix. Not modified.
     * @param upper true of upper triangular and false for lower.
     * @param tol How close to being identical each element needs to be.
     * @return true if equals and false otherwise.
     */
    public static boolean isEqualsTriangle(ReshapeMatrix32F a, ReshapeMatrix32F b, boolean upper, float tol)
    {
        if( a.numRows != b.numRows || a.numCols != b.numCols ) {
            return false;
        }

        if( upper ) {
            for( int i = 0; i < a.numRows; i++ ) {
                for( int j = i; j < a.numCols; j++ ) {
                    if( Math.abs(a.get(i,j)-b.get(i,j)) > tol )
                        return false;
                }
            }
        } else {
            for( int i = 0; i < a.numRows; i++ ) {
                int end = Math.min(i,a.numCols-1);

                for( int j = 0; j <= end; j++ ) {
                    if( Math.abs(a.get(i,j)-b.get(i,j)) > tol )
                        return false;
                }
            }
        }

        return true;
    }

    /**
     * <p>
     * Checks to see if each element in the two matrices are equal:
     * a<sub>ij</sub> == b<sub>ij</sub>
     * <p>
     *
     * <p>
     * NOTE: If any of the elements are NaN then false is returned.  If two corresponding
     * elements are both positive or negative infinity then they are equal.
     * </p>
     * 
     * @param a A matrix. Not modified.
     * @param b A matrix. Not modified.
     * @return true if identical and false otherwise.
     */
    public static boolean isEquals( D1Matrix32F a, D1Matrix32F b ) {
        if( a.numRows != b.numRows || a.numCols != b.numCols ) {
            return false;
        }

        final int length = a.getNumElements();
        for( int i = 0; i < length; i++ ) {
            if( !(a.get(i) == b.get(i)) ) {
                return false;
            }
        }

        return true;
    }

    /**
     * <p>
     * Checks to see if each corresponding element in the two matrices are
     * within tolerance of each other or have the some symbolic meaning.  This
     * can handle NaN and Infinite numbers.
     * <p>
     *
     * <p>
     * If both elements are countable then the following equality test is used:<br>
     * |a<sub>ij</sub> - b<sub>ij</sub>| &le; tol.<br>
     * Otherwise both numbers must both be float.NaN, float.POSITIVE_INFINITY, or
     * float.NEGATIVE_INFINITY to be identical.
     * </p>
     *
     * @param a A matrix. Not modified.
     * @param b A matrix. Not modified.
     * @param tol Tolerance for equality.
     * @return true if identical and false otherwise.
     */
     public static boolean isIdentical( D1Matrix32F a, D1Matrix32F b , float tol ) {
        if( a.numRows != b.numRows || a.numCols != b.numCols ) {
            return false;
        }
        if( tol < 0 )
            throw new IllegalArgumentException("Tolerance must be greater than or equal to zero.");

        final int length = a.getNumElements();
        for( int i = 0; i < length; i++ ) {
            float valA = a.get(i);
            float valB = b.get(i);

            // if either is negative or positive infinity the result will be positive infinity
            // if either is NaN the result will be NaN
            float diff = Math.abs(valA-valB);

            // diff = NaN == false
            // diff = infinity == false
            if( tol >= diff )
                continue;

            if( Float.isNaN(valA) ) {
                return Float.isNaN(valB);
            } else if( Float.isInfinite(valA) ) {
                return valA == valB;
            } else {
                return false;
            }
        }

        return true;
    }

    /**
     * <p>
     * Checks to see if a matrix is orthogonal or isometric.
     * </p>
     *
     * @param Q The matrix being tested. Not modified.
     * @param tol Tolerance.
     * @return True if it passes the test.
     */
    public static boolean isOrthogonal( DenseMatrix32F Q , float tol )
    {
       if( Q.numRows < Q.numCols ) {
            throw new IllegalArgumentException("The number of rows must be more than or equal to the number of columns");
        }

        DenseMatrix32F u[] = CommonOps.columnsToVector(Q, null);

        for( int i = 0; i < u.length; i++ ) {
            DenseMatrix32F a = u[i];

            for( int j = i+1; j < u.length; j++ ) {
                float val = VectorVectorMult.innerProd(a,u[j]);

                if( !(Math.abs(val) <= tol))
                    return false;
            }
        }

        return true;
    }

    /**
     * Checks to see if the rows of the provided matrix are linearly independent.
     *
     * @param A Matrix whose rows are being tested for linear independence.
     * @return true if linearly independent and false otherwise.
     */
    public static boolean isRowsLinearIndependent( DenseMatrix32F A )
    {
        // LU decomposition
        LUDecomposition<DenseMatrix32F> lu = DecompositionFactory.lu(A.numRows,A.numCols);
        if( lu.inputModified() )
            A = A.copy();

        if( !lu.decompose(A))
            throw new RuntimeException("Decompositon failed?");

        // if they are linearly independent it should not be singular
        return !lu.isSingular();
    }

    /**
     * Checks to see if the provided matrix is within tolerance to an identity matrix.
     *
     * @param mat Matrix being examined.  Not modified.
     * @param tol Tolerance.
     * @return True if it is within tolerance to an identify matrix.
     */
    public static boolean isIdentity( DenseMatrix32F mat , float tol )
    {
        // see if the result is an identity matrix
        int index = 0;
        for( int i = 0; i < mat.numRows; i++ ) {
            for( int j = 0; j < mat.numCols; j++ ) {
                if( i == j ) {
                    if( !(Math.abs(mat.get(index++)-1) <= tol) )
                        return false;
                } else {
                    if( !(Math.abs(mat.get(index++)) <= tol) )
                        return false;
                }
            }
        }

        return true;
    }

    /**
     * Checks to see if every value in the matrix is the specified value.
     *
     * @param mat The matrix being tested.  Not modified.
     * @param val Checks to see if every element in the matrix has this value.
     * @param tol True if all the elements are within this tolerance.
     * @return true if the test passes.
     */
    public static boolean isConstantVal( DenseMatrix32F mat , float val , float tol )
    {
        // see if the result is an identity matrix
        int index = 0;
        for( int i = 0; i < mat.numRows; i++ ) {
            for( int j = 0; j < mat.numCols; j++ ) {
                if( !(Math.abs(mat.get(index++)-val) <= tol) )
                    return false;

            }
        }

        return true;
    }

    /**
     * Checks to see if all the diagonal elements in the matrix are positive.
     *
     * @param a A matrix. Not modified.
     * @return true if all the  diagonal elements are positive, false otherwise.
     */
    public static boolean isDiagonalPositive( DenseMatrix32F a ) {
        for( int i = 0; i < a.numRows; i++ ) {
            if( !(a.get(i,i) >= 0) )
                return false;
        }
        return true;
    }

    // TODO write this
    public static boolean isFullRank( DenseMatrix32F a ) {
        throw new RuntimeException("Implement");
    }

    /**
     * <p>
     * Checks to see if the two matrices are the negative of each other:<br>
     * <br>
     * a<sub>ij</sub> = -b<sub>ij</sub>
     * </p>
     *
     * @param a First matrix.  Not modified.
     * @param b Second matrix.  Not modified.
     * @param tol Numerical tolerance.
     * @return True if they are the negative of each other within tolerance.
     */
    public static boolean isNegative(D1Matrix32F a, D1Matrix32F b, float tol) {
        if( a.numRows != b.numRows || a.numCols != b.numCols )
            throw new IllegalArgumentException("Matrix dimensions must match");

        int length = a.getNumElements();

        for( int i = 0; i < length; i++ ) {
            if( !(Math.abs(a.get(i)+b.get(i)) <= tol) )
                return false;
        }

        return true;
    }

    /**
     * <p>
     * Checks to see if a matrix is upper triangular or Hessenberg. A Hessenberg matrix of degree N
     * has the following property:<br>
     * <br>
     * a<sub>ij</sub> &le; 0 for all i < j+N<br>
     * <br>
     * A triangular matrix is a Hessenberg matrix of degree 0.
     * </p>
     * @param A Matrix being tested.  Not modified.
     * @param hessenberg The degree of being hessenberg.
     * @param tol How close to zero the lower left elements need to be.
     * @return If it is an upper triangular/hessenberg matrix or not.
     */
    public static boolean isUpperTriangle(DenseMatrix32F A , int hessenberg , float tol ) {
        if( A.numRows != A.numCols )
            return false;

        for( int i = hessenberg+1; i < A.numRows; i++ ) {
            for( int j = 0; j < i-hessenberg; j++ ) {
                if( !(Math.abs(A.get(i,j)) <= tol) ) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Computes the rank of a matrix using a default tolerance.
     *
     * @param A Matrix whose rank is to be calculated.  Not modified.
     * @return The matrix's rank.
     */
    public static int rank( DenseMatrix32F A ) {
        return rank(A, UtilEjml.EPS*100);
    }

    /**
     * Computes the rank of a matrix using the specified tolerance.
     *
     * @param A Matrix whose rank is to be calculated.  Not modified.
     * @param threshold The numerical threshold used to determine a singular value.
     * @return The matrix's rank.
     */
    public static int rank( DenseMatrix32F A , float threshold ) {
        SingularValueDecomposition<DenseMatrix32F> svd = DecompositionFactory.svd(A.numRows,A.numCols,false,false,true);

        if( svd.inputModified() )
            A = A.copy();

        if( !svd.decompose(A) )
            throw new RuntimeException("Decomposition failed");

        return SingularOps.rank(svd, threshold);
    }

    /**
     * Computes the nullity of a matrix using the default tolerance. 
     *
     * @param A Matrix whose rank is to be calculated.  Not modified.
     * @return The matrix's nullity.
     */
    public static int nullity( DenseMatrix32F A ) {
        return nullity(A, UtilEjml.EPS*100);
    }

    /**
     * Computes the nullity of a matrix using the specified tolerance.
     *
     * @param A Matrix whose rank is to be calculated.  Not modified.
     * @param threshold The numerical threshold used to determine a singular value.
     * @return The matrix's nullity.
     */
    public static int nullity( DenseMatrix32F A , float threshold ) {
        SingularValueDecomposition<DenseMatrix32F> svd = DecompositionFactory.svd(A.numRows,A.numCols,false,false,true);

        if( svd.inputModified() )
            A = A.copy();

        if( !svd.decompose(A) )
            throw new RuntimeException("Decomposition failed");

        return SingularOps.nullity(svd,threshold);
    }
}
