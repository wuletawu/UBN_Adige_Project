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
package org.jgrasstools.hortonmachine.models.hm;

import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.hortonmachine.modules.basin.topindex.OmsTopIndex;
import org.jgrasstools.hortonmachine.utils.HMTestCase;
import org.jgrasstools.hortonmachine.utils.HMTestMaps;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
/**
 * Test {@link OmsTopIndex}.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class TestTopIndex extends HMTestCase {

    public void testAb() throws Exception {
        HashMap<String, Double> envelopeParams = HMTestMaps.getEnvelopeparams();
        CoordinateReferenceSystem crs = HMTestMaps.getCrs();

        double[][] tcaData = HMTestMaps.tcaData;
        GridCoverage2D tcaCoverage = CoverageUtilities.buildCoverage("tca", tcaData, envelopeParams, crs, true);
        double[][] slopeData = HMTestMaps.slopeData;
        GridCoverage2D slopeCoverage = CoverageUtilities.buildCoverage("plan", slopeData, envelopeParams, crs, true);

        OmsTopIndex topindex = new OmsTopIndex();
        topindex.inTca = tcaCoverage;
        topindex.inSlope = slopeCoverage;
        topindex.pm = pm;
        topindex.process();

        GridCoverage2D topindexCoverage = topindex.outTopindex;

        checkMatrixEqual(topindexCoverage.getRenderedImage(), HMTestMaps.topIndexData, 0.01);
    }

}
