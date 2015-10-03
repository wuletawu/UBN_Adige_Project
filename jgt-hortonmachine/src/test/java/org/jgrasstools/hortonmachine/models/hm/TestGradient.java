package org.jgrasstools.hortonmachine.models.hm;

import java.io.IOException;
import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.hortonmachine.modules.geomorphology.gradient.OmsGradient;
import org.jgrasstools.hortonmachine.utils.HMTestCase;
import org.jgrasstools.hortonmachine.utils.HMTestMaps;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
/**
 * It test the {@link OmsGradient} module.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class TestGradient extends HMTestCase {

    public void testGradient() throws IOException {

        HashMap<String, Double> envelopeParams = HMTestMaps.getEnvelopeparams();
        CoordinateReferenceSystem crs = HMTestMaps.getCrs();
        double[][] pitData = HMTestMaps.pitData;
        GridCoverage2D pitfillerCoverage = CoverageUtilities.buildCoverage("elevation", pitData, envelopeParams, crs, true);

        OmsGradient gradient = new OmsGradient();
        gradient.inElev = pitfillerCoverage;
        gradient.pm = pm;

        gradient.process();

        GridCoverage2D gradientCoverage = gradient.outSlope;
        checkMatrixEqual(gradientCoverage.getRenderedImage(), HMTestMaps.gradientData, 0.01);
    }

    public void testGradientHorn() throws IOException {

        HashMap<String, Double> envelopeParams = HMTestMaps.getEnvelopeparams();
        CoordinateReferenceSystem crs = HMTestMaps.getCrs();
        double[][] pitData = HMTestMaps.pitData;
        GridCoverage2D pitfillerCoverage = CoverageUtilities.buildCoverage("elevation", pitData, envelopeParams, crs, true);

        OmsGradient gradient = new OmsGradient();
        gradient.inElev = pitfillerCoverage;
        gradient.pm = pm;
        gradient.pMode = 1;

        gradient.process();

        GridCoverage2D gradientCoverage = gradient.outSlope;
        checkMatrixEqual(gradientCoverage.getRenderedImage(), HMTestMaps.gradientHornData, 0.01);
    }

    public void testGradientEvans() throws IOException {

        HashMap<String, Double> envelopeParams = HMTestMaps.getEnvelopeparams();
        CoordinateReferenceSystem crs = HMTestMaps.getCrs();
        double[][] pitData = HMTestMaps.pitData;
        GridCoverage2D pitfillerCoverage = CoverageUtilities.buildCoverage("elevation", pitData, envelopeParams, crs, true);

        OmsGradient gradient = new OmsGradient();
        gradient.inElev = pitfillerCoverage;
        gradient.pm = pm;
        gradient.pMode = 2;

        gradient.process();

        GridCoverage2D gradientCoverage = gradient.outSlope;
        checkMatrixEqual(gradientCoverage.getRenderedImage(), HMTestMaps.gradientEvansData, 0.01);
    }

}
