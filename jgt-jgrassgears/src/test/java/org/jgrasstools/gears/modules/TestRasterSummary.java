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
package org.jgrasstools.gears.modules;

import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.modules.r.summary.OmsRasterSummary;
import org.jgrasstools.gears.utils.HMTestCase;
import org.jgrasstools.gears.utils.HMTestMaps;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

/**
 * Test for {@link OmsRasterSummary}.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class TestRasterSummary extends HMTestCase {
    public void testCoverageSummary() throws Exception {

        double[][] inData = HMTestMaps.extractNet0Data;
        HashMap<String, Double> envelopeParams = HMTestMaps.getEnvelopeparams();
        CoordinateReferenceSystem crs = HMTestMaps.getCrs();
        GridCoverage2D inCoverage = CoverageUtilities.buildCoverage("data", inData, envelopeParams, crs, true);

        OmsRasterSummary summary = new OmsRasterSummary();
        summary.pm = pm;
        summary.inRaster = inCoverage;
        summary.doHistogram = true;
        summary.pBins = 100;
        summary.process();

        double min = summary.outMin;
        double max = summary.outMax;
        double mean = summary.outMean;
        double sdev = summary.outSdev;
        double range = summary.outRange;
        double sum = summary.outSum;

        assertEquals(2.0, min);
        assertEquals(2.0, max);
        assertEquals(2.0, mean);
        assertEquals(0.0, sdev);
        assertEquals(0.0, range);
        assertEquals(18.0, sum);

        double[] minMax = OmsRasterSummary.getMinMax(inCoverage);
        assertEquals(2.0, minMax[0]);
        assertEquals(2.0, minMax[1]);

        double[][] cb = summary.outCb;
        // for( int i = 0; i < cb.length; i++ ) {
        // System.out.println(cb[i][0] + "\t" + cb[i][1] + "\t" + cb[i][2] + "%");
        // }

        assertEquals(cb[0][0], 2.0);
        assertEquals(cb[0][1], 9.0);
        assertEquals(cb[0][2], 11.25);
        assertTrue(JGTConstants.isNovalue(cb[cb.length - 1][0]));
        assertEquals(cb[cb.length - 1][1], 71.00);
        assertEquals(cb[cb.length - 1][2], 88.75);

        // System.out.println();

    }

}
