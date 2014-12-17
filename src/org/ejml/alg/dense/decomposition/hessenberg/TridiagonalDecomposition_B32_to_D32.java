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

package org.ejml.alg.dense.decomposition.hessenberg;

import org.ejml.EjmlParameters;
import org.ejml.alg.block.decomposition.hessenberg.TridiagonalDecompositionHouseholder_B32;
import org.ejml.alg.dense.decomposition.BaseDecomposition_B32_to_D32;
import org.ejml.data.BlockMatrix32F;
import org.ejml.data.DenseMatrix32F;
import org.ejml.interfaces.decomposition.TridiagonalSimilarDecomposition;
import org.ejml.ops.CommonOps;


/**
 * Wrapper around a block implementation of TridiagonalSimilarDecomposition
 *
 * @author Peter Abeles
 */
public class TridiagonalDecomposition_B32_to_D32
        extends BaseDecomposition_B32_to_D32
        implements TridiagonalSimilarDecomposition<DenseMatrix32F> {


    public TridiagonalDecomposition_B32_to_D32() {
        this(EjmlParameters.BLOCK_WIDTH);
    }

    public TridiagonalDecomposition_B32_to_D32(int blockSize) {
        super(new TridiagonalDecompositionHouseholder_B32(),blockSize);
    }

    @Override
    public DenseMatrix32F getT(DenseMatrix32F T) {
        int N = Ablock.numRows;

        if( T == null ) {
            T = new DenseMatrix32F(N,N);
        } else {
            CommonOps.fill(T, 0);
        }

        float[] diag = new float[ N ];
        float[] off = new float[ N ];

        ((TridiagonalDecompositionHouseholder_B32)alg).getDiagonal(diag,off);

        T.unsafe_set(0,0,diag[0]);
        for( int i = 1; i < N; i++ ) {
            T.unsafe_set(i,i,diag[i]);
            T.unsafe_set(i,i-1,off[i-1]);
            T.unsafe_set(i-1,i,off[i-1]);
        }

        return T;
    }

    @Override
    public DenseMatrix32F getQ(DenseMatrix32F Q, boolean transposed) {
        if( Q == null ) {
            Q = new DenseMatrix32F(Ablock.numRows,Ablock.numCols);
        }

        BlockMatrix32F Qblock = new BlockMatrix32F();
        Qblock.numRows =  Q.numRows;
        Qblock.numCols =  Q.numCols;
        Qblock.blockLength = blockLength;
        Qblock.data = Q.data;

        ((TridiagonalDecompositionHouseholder_B32)alg).getQ(Qblock,transposed);

        convertBlockToRow(Q.numRows,Q.numCols,Ablock.blockLength,Q.data);

        return Q;
    }

    @Override
    public void getDiagonal(float[] diag, float[] off) {
        ((TridiagonalDecompositionHouseholder_B32)alg).getDiagonal(diag,off);
    }
}
