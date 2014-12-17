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

package org.ejml.alg.dense.decomposition;

import org.ejml.alg.block.BlockMatrixOps;
import org.ejml.data.BlockMatrix32F;
import org.ejml.data.DenseMatrix32F;
import org.ejml.interfaces.decomposition.DecompositionInterface;


/**
 * Generic interface for wrapping a {@link BlockMatrix32F} decomposition for
 * processing of {@link org.ejml.data.DenseMatrix32F}.
 *
 * @author Peter Abeles
 */
public class BaseDecomposition_B32_to_D32 implements DecompositionInterface<DenseMatrix32F> {

    protected DecompositionInterface<BlockMatrix32F> alg;

    protected float[]tmp;
    protected BlockMatrix32F Ablock = new BlockMatrix32F();
    protected int blockLength;

    public BaseDecomposition_B32_to_D32(DecompositionInterface<BlockMatrix32F> alg,
                                        int blockLength) {
        this.alg = alg;
        this.blockLength = blockLength;
    }

    @Override
    public boolean decompose(DenseMatrix32F A) {
        Ablock.numRows = A.numRows;
        Ablock.numCols = A.numCols;
        Ablock.blockLength = blockLength;
        Ablock.data = A.data;

        int tmpLength = Math.min( Ablock.blockLength , A.numRows ) * A.numCols;

        if( tmp == null || tmp.length < tmpLength )
            tmp = new float[ tmpLength ];

        // doing an in-place convert is much more memory efficient at the cost of a little
        // but of CPU
        BlockMatrixOps.convertRowToBlock(A.numRows,A.numCols,Ablock.blockLength,A.data,tmp);

        boolean ret = alg.decompose(Ablock);

        // convert it back to the normal format if it wouldn't have been modified
        if( !alg.inputModified() ) {
            BlockMatrixOps.convertBlockToRow(A.numRows,A.numCols,Ablock.blockLength,A.data,tmp);
        }

        return ret;
    }

    public void convertBlockToRow(int numRows , int numCols , int blockLength ,
                                  float[] data) {
        int tmpLength = Math.min( blockLength , numRows ) * numCols;

        if( tmp == null || tmp.length < tmpLength )
            tmp = new float[ tmpLength ];

        BlockMatrixOps.convertBlockToRow(numRows,numCols,Ablock.blockLength,data,tmp);
    }

    @Override
    public boolean inputModified() {
        return alg.inputModified();
    }
}
