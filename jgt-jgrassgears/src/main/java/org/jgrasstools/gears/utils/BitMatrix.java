/*
 * This file is part of JGrasstools (http://www.jgrasstools.org)
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * JGrasstools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.jgrasstools.gears.utils;

import java.util.BitSet;

/**
 * A matrix object that is backed by a bitmap.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class BitMatrix {

    private BitSet bitSet;
    private int cols;

    public BitMatrix( int cols, int rows ) {
        this.cols = cols;
        bitSet = new BitSet(cols * rows);
    }

    /**
     * Gets the state of a particular bit.
     * 
     * @param col the col in the matrix of the bit.
     * @param row the row in the matrix of the bit.
     * @return <code>true</code> if the bit is set, else <code>false</code>.
     */
    public boolean isMarked( int col, int row ) {
        return bitSet.get(row * cols + col);
    }

    /**
     * Marks the bit in a given position.
     * 
     * @param col the col position to mark.
     * @param row the row position to mark.
     */
    public void mark( int col, int row ) {
        bitSet.set(row * cols + col);
    }

    /**
     * Unmarks the bit in a given position.
     * 
     * @param col the col position to mark.
     * @param row the row position to mark.
     */
    public void unMark( int col, int row ) {
        bitSet.set(row * cols + col, false);
    }
}
