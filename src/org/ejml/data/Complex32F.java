/*
 * Copyright (c) 2009-2013, Peter Abeles. All Rights Reserved.
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

package org.ejml.data;

import java.io.Serializable;

/**
 * <p>
 * Represents a complex number using 64bit floating point numbers.  A complex number is composed of
 * a real and imaginary components.
 * </p>
 */
public class Complex32F implements Serializable {
    public float real;
    public float imaginary;

    public Complex32F(float real, float imaginary) {
        this.real = real;
        this.imaginary = imaginary;
    }

    public Complex32F() {
    }

    public float getReal() {
        return real;
    }

    public float getMagnitude() {
        return (float)Math.sqrt(real*real + imaginary*imaginary);
    }

    public float getMagnitude2() {
        return real*real + imaginary*imaginary;
    }

    public void setReal(float real) {
        this.real = real;
    }

    public float getImaginary() {
        return imaginary;
    }

    public void setImaginary(float imaginary) {
        this.imaginary = imaginary;
    }

    public void set(float real, float imaginary) {
        this.real = real;
        this.imaginary = imaginary;
    }

    public boolean isReal() {
        return imaginary == 0.0;
    }

    public String toString() {
        if( imaginary == 0 ) {
            return ""+real;
        } else {
            return real+" "+imaginary+"i";
        }
    }
}
