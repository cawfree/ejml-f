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

package org.ejml.alg.dense.decomposition.chol;

import org.ejml.EjmlParameters;
import org.ejml.alg.block.BlockMatrixOps;
import org.ejml.alg.block.decomposition.chol.CholeskyOuterForm_B32;
import org.ejml.alg.dense.decomposition.BaseDecomposition_B32_to_D32;
import org.ejml.data.BlockMatrix32F;
import org.ejml.data.DenseMatrix32F;
import org.ejml.interfaces.decomposition.CholeskyDecomposition;


/**
 * Wrapper around {@link org.ejml.alg.block.decomposition.chol.CholeskyOuterForm_B32} that allows
 * it to process DenseMatrix32F.
 *
 * @author Peter Abeles
 */
public class CholeskyDecomposition_B32_to_D32
        extends BaseDecomposition_B32_to_D32 implements CholeskyDecomposition<DenseMatrix32F> {

    public CholeskyDecomposition_B32_to_D32(boolean lower) {
        super(new CholeskyOuterForm_B32(lower), EjmlParameters.BLOCK_WIDTH);
    }

    @Override
    public boolean isLower() {
        return ((CholeskyOuterForm_B32)alg).isLower();
    }

    @Override
    public DenseMatrix32F getT(DenseMatrix32F T) {
        BlockMatrix32F T_block = ((CholeskyOuterForm_B32)alg).getT(null);

        if( T == null ) {
            T = new DenseMatrix32F(T_block.numRows,T_block.numCols);
        }

        BlockMatrixOps.convert(T_block,T);
        // todo set zeros
        return T;
    }
}
