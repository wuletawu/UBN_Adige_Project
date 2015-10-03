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
package org.jgrasstools.gears.modules.r.summary;

import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_DOCUMENTATION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_doHistogram_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_inRaster_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_outCb_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_outMax_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_outMean_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_outMin_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_outRange_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_outSdev_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_outSum_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERSUMMARY_pBins_DESCRIPTION;

import java.awt.image.RenderedImage;
import java.util.List;

import javax.media.jai.JAI;
import javax.media.jai.ParameterBlockJAI;
import javax.media.jai.RenderedOp;

import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jaitools.media.jai.zonalstats.Result;
import org.jaitools.media.jai.zonalstats.ZonalStats;
import org.jaitools.media.jai.zonalstats.ZonalStatsDescriptor;
import org.jaitools.numeric.Statistic;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.utils.math.CoupledFieldsMoments;

@Description(OMSRASTERSUMMARY_DESCRIPTION)
@Documentation(OMSRASTERSUMMARY_DOCUMENTATION)
@Author(name = OMSRASTERSUMMARY_AUTHORNAMES, contact = OMSRASTERSUMMARY_AUTHORCONTACTS)
@Keywords(OMSRASTERSUMMARY_KEYWORDS)
@Label(OMSRASTERSUMMARY_LABEL)
@Name(OMSRASTERSUMMARY_NAME)
@Status(OMSRASTERSUMMARY_STATUS)
@License(OMSRASTERSUMMARY_LICENSE)
public class OmsRasterSummary extends JGTModel {

    @Description(OMSRASTERSUMMARY_inRaster_DESCRIPTION)
    @In
    public GridCoverage2D inRaster;

    @Description(OMSRASTERSUMMARY_pBins_DESCRIPTION)
    @In
    public int pBins = 100;

    @Description(OMSRASTERSUMMARY_doHistogram_DESCRIPTION)
    @In
    public boolean doHistogram = false;

    @Description(OMSRASTERSUMMARY_outMin_DESCRIPTION)
    @Out
    public Double outMin = null;

    @Description(OMSRASTERSUMMARY_outMax_DESCRIPTION)
    @Out
    public Double outMax = null;

    @Description(OMSRASTERSUMMARY_outMean_DESCRIPTION)
    @Out
    public Double outMean = null;

    @Description(OMSRASTERSUMMARY_outSdev_DESCRIPTION)
    @Out
    public Double outSdev = null;

    @Description(OMSRASTERSUMMARY_outRange_DESCRIPTION)
    @Out
    public Double outRange = null;

    @Description(OMSRASTERSUMMARY_outSum_DESCRIPTION)
    @Out
    public Double outSum = null;

    @Description(OMSRASTERSUMMARY_outCb_DESCRIPTION)
    @Out
    public double[][] outCb = null;

    private Statistic[] stats;

    @Execute
    public void process() throws Exception {
        if (!concatOr(outMin == null, doReset)) {
            return;
        }

        // TODO use the geotools bridge instead of jaitools:
        // http://svn.osgeo.org/geotools/trunk/modules/library/coverage/src/test/java/org/geotools/coverage/processing/operation/ZonalStasTest.java

        RenderedImage inRI = inRaster.getRenderedImage();
        ParameterBlockJAI pb = new ParameterBlockJAI("ZonalStats");
        pb.setSource("dataImage", inRI);
        // pb.setSource("zoneImage", null);

        if (stats == null) {
            stats = new Statistic[]{Statistic.MIN, Statistic.MAX, Statistic.MEAN, Statistic.SDEV, Statistic.RANGE, Statistic.SUM};
        }
        pb.setParameter("stats", stats);

        RenderedOp op = JAI.create("ZonalStats", pb);

        ZonalStats zonalStats = (ZonalStats) op.getProperty(ZonalStatsDescriptor.ZONAL_STATS_PROPERTY);
        List<Result> results = zonalStats.results();
        for( Result result : results ) {
            Statistic statistic = result.getStatistic();
            Double value = result.getValue();

            switch( statistic ) {
            case MIN:
                outMin = value;
                break;
            case MAX:
                outMax = value;
                break;
            case MEAN:
                outMean = value;
                break;
            case SDEV:
                outSdev = value;
                break;
            case RANGE:
                outRange = value;
                break;
            case SUM:
                outSum = value;
                break;
            default:
                break;
            }
        }

        if (!doHistogram)
            return;

        double[][] cb = new CoupledFieldsMoments().process(inRI, null, pBins, 1, 2, pm, 1);

        int width = inRI.getWidth();
        int height = inRI.getHeight();
        int pixelsNum = width * height;
        outCb = new double[cb.length + 1][3];

        double sum = 0;
        for( int i = 0; i < outCb.length; i++ ) {
            if (i < outCb.length - 1) {
                outCb[i][0] = cb[i][0];
                outCb[i][1] = cb[i][1];
                sum = sum + cb[i][1];
                outCb[i][2] = cb[i][1] * 100.0 / pixelsNum;
            } else {
                outCb[i][0] = JGTConstants.doubleNovalue;
                double nans = pixelsNum - sum;
                outCb[i][1] = nans;
                outCb[i][2] = nans * 100.0 / pixelsNum;
            }

        }

    }

    public static double[] getMinMax( GridCoverage2D raster ) throws Exception {
        OmsRasterSummary summary = new OmsRasterSummary();
        summary.inRaster = raster;
        summary.doHistogram = false;
        summary.stats = new Statistic[]{Statistic.MIN, Statistic.MAX};
        summary.process();

        double min = summary.outMin;
        double max = summary.outMax;
        return new double[]{min, max};
    }

    public static double[] getMinMaxAvgSum( GridCoverage2D raster ) throws Exception {
        OmsRasterSummary summary = new OmsRasterSummary();
        summary.inRaster = raster;
        summary.doHistogram = false;
        summary.stats = new Statistic[]{Statistic.MIN, Statistic.MAX, Statistic.MEAN, Statistic.SUM};
        summary.process();

        double min = summary.outMin;
        double max = summary.outMax;
        double avg = summary.outMean;
        double sum = summary.outSum;
        return new double[]{min, max, sum, avg};
    }

}
