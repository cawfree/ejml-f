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

package org.ejml.alg.dense.decomposition.eig;

import org.ejml.alg.dense.decomposition.eig.watched.WatchedTwoStepQREigenvalue;
import org.ejml.alg.dense.decomposition.eig.watched.WatchedTwoStepQREigenvector;
import org.ejml.alg.dense.decomposition.hessenberg.HessenbergSimilarDecomposition_D32;
import org.ejml.data.Complex32F;
import org.ejml.data.DenseMatrix32F;
import org.ejml.interfaces.decomposition.EigenDecomposition;


/**
 * <p>
 * Finds the eigenvalue decomposition of an arbitrary square matrix using the implicit float-step QR algorithm.
 * Watched is included in its name because it is designed to print out internal debugging information.  This
 * class is still underdevelopment and has yet to be optimized.
 * </p>
 *
 * <p>
 * Based off the description found in:<br>
 * David S. Watkins, "Fundamentals of Matrix Computations." Second Edition.
 * </p>
 *
 * @author Peter Abeles
 */
//TODO looks like there might be some pointless copying of arrays going on
public class WatchedTwoStepQRDecomposition_D32
        implements EigenDecomposition<DenseMatrix32F> {

    HessenbergSimilarDecomposition_D32 hessenberg;
    WatchedTwoStepQREigenvalue algValue;
    WatchedTwoStepQREigenvector algVector;

    DenseMatrix32F H;

    // should it compute eigenvectors or just eigenvalues
    boolean computeVectors;

    public WatchedTwoStepQRDecomposition_D32(boolean computeVectors) {
        hessenberg = new HessenbergSimilarDecomposition_D32(10);
        algValue = new WatchedTwoStepQREigenvalue();
        algVector = new WatchedTwoStepQREigenvector();

        this.computeVectors = computeVectors;
    }

    @Override
    public boolean decompose(DenseMatrix32F A) {

        if( !hessenberg.decompose(A) )
            return false;

        H = hessenberg.getH(null);

        algValue.getImplicitQR().createR = false;
//        algValue.getImplicitQR().setChecks(true,true,true);

        if( !algValue.process(H) )
            return false;

//        for( int i = 0; i < A.numRows; i++ ) {
//            System.out.println(algValue.getEigenvalues()[i]);
//        }

        algValue.getImplicitQR().createR = true;

        if( computeVectors )
            return algVector.process(algValue.getImplicitQR(), H, hessenberg.getQ(null));
        else
            return true;
    }

    @Override
    public boolean inputModified() {
        return hessenberg.inputModified();
    }

    @Override
    public int getNumberOfEigenvalues() {
        return algValue.getEigenvalues().length;
    }

    @Override
    public Complex32F getEigenvalue(int index) {
        return algValue.getEigenvalues()[index];
    }

    @Override
    public DenseMatrix32F getEigenVector(int index) {
        return algVector.getEigenvectors()[index];
    }
}
