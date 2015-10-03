/*
 * JGrass - Free Open Source Java GIS http://www.jgrass.org 
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Library General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Library General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU Library General Public License
 * along with this library; if not, write to the Free Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
package org.jgrasstools.gears.utils.math.interpolation;

import java.util.List;

import org.jgrasstools.gears.libs.modules.JGTConstants;

/**
 * A polynomial interpolation function.
 * 
 * <p>
 * This was done basing on the Java Number Cruncher Book by Ronald Mak
 * (see http://www.apropos-logic.com/books.html).
 * </p>
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class PolynomialInterpolator implements Interpolator {
    private int n;
    private double[][] data;
    /** 
     * divided difference table 
     */
    private double dd[][];

    /**
     * Constructor.
     * 
     * @param xList the list of X samples.
     * @param yList the list of Y = f(X) samples.
     */
    public PolynomialInterpolator( List<Double> xList, List<Double> yList ) {
        this.dd = new double[xList.size()][xList.size()];

        data = new double[xList.size()][xList.size()];
        for( int i = 0; i < xList.size(); i++ ) {
            data[i][0] = xList.get(i);
            data[i][1] = yList.get(i);
        }

        for( int i = 0; i < data.length; ++i ) {
            addPoint(data[i]);
        }
    }

    private void addPoint( double[] dataPoint ) {
        if (n >= data.length)
            return;

        // data[n] = dataPoint;
        dd[n][0] = dataPoint[1];

        ++n;

        for( int order = 1; order < n; ++order ) {
            int bottom = n - order - 1;
            double numerator = dd[bottom + 1][order - 1] - dd[bottom][order - 1];
            double denominator = data[bottom + order][0] - data[bottom][0];

            dd[bottom][order] = numerator / denominator;
        }
    }

    public double getInterpolated( double x ) {
        if (n < 2)
            return JGTConstants.doubleNovalue;

        double y = dd[0][0];
        double xFactor = 1;

        // Compute the value of the function.
        for( int order = 1; order < n; ++order ) {
            xFactor = xFactor * (x - data[order - 1][0]);
            y = y + xFactor * dd[0][order];
        }
        return y;
    }

}