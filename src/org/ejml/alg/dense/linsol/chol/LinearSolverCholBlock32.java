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

package org.ejml.alg.dense.linsol.chol;

import org.ejml.alg.block.BlockMatrixOps;
import org.ejml.alg.block.linsol.chol.BlockCholeskyOuterSolver;
import org.ejml.alg.dense.linsol.LinearSolver_B32_to_D32;
import org.ejml.data.DenseMatrix32F;


/**
 * A wrapper around {@link org.ejml.interfaces.decomposition.CholeskyDecomposition}(BlockMatrix32F) that allows
 * it to be easily used with {@link org.ejml.data.DenseMatrix32F}.
 *
 * @author Peter Abeles
 */
public class LinearSolverCholBlock32 extends LinearSolver_B32_to_D32 {

    public LinearSolverCholBlock32() {
        super(new BlockCholeskyOuterSolver());
    }

    /**
     * Only converts the B matrix and passes that onto solve.  Te result is then copied into
     * the input 'X' matrix.
     * 
     * @param B A matrix &real; <sup>m &times; p</sup>.  Not modified.
     * @param X A matrix &real; <sup>n &times; p</sup>, where the solution is written to.  Modified.
     */
    @Override
    public void solve(DenseMatrix32F B, DenseMatrix32F X) {
        blockB.reshape(B.numRows,B.numCols,false);
        BlockMatrixOps.convert(B,blockB);

        // since overwrite B is true X does not need to be passed in
        alg.solve(blockB,null);

        BlockMatrixOps.convert(blockB,X);
    }

}
