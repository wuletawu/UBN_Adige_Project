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
package org.jgrasstools.gears.modules.r.cutout;

import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_DOCUMENTATION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_doInverse_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_inMask_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_inRaster_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_outRaster_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_pMax_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCUTOUT_pMin_DESCRIPTION;
import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;

import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
import java.util.HashMap;

import javax.media.jai.iterator.RandomIter;
import javax.media.jai.iterator.RandomIterFactory;
import javax.media.jai.iterator.WritableRandomIter;

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
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;

@Description(OMSCUTOUT_DESCRIPTION)
@Documentation(OMSCUTOUT_DOCUMENTATION)
@Author(name = OMSCUTOUT_AUTHORNAMES, contact = OMSCUTOUT_AUTHORCONTACTS)
@Keywords(OMSCUTOUT_KEYWORDS)
@Label(OMSCUTOUT_LABEL)
@Name(OMSCUTOUT_NAME)
@Status(OMSCUTOUT_STATUS)
@License(OMSCUTOUT_LICENSE)
public class OmsCutOut extends JGTModel {

    @Description(OMSCUTOUT_inRaster_DESCRIPTION)
    @In
    public GridCoverage2D inRaster;

    @Description(OMSCUTOUT_inMask_DESCRIPTION)
    @In
    public GridCoverage2D inMask;

    @Description(OMSCUTOUT_pMax_DESCRIPTION)
    @In
    public Double pMax;

    @Description(OMSCUTOUT_pMin_DESCRIPTION)
    @In
    public Double pMin;

    @Description(OMSCUTOUT_doInverse_DESCRIPTION)
    @In
    public boolean doInverse = false;

    @Description(OMSCUTOUT_outRaster_DESCRIPTION)
    @Out
    public GridCoverage2D outRaster = null;

    private RandomIter maskIter;
    private boolean doMax = false;
    private boolean doMin = false;

    @Execute
    public void process() throws Exception {
        if (!concatOr(outRaster == null, doReset)) {
            return;
        }

        double max = -1;
        double min = -1;
        // do autoboxing only once
        if (pMax != null) {
            max = pMax;
            doMax = true;
        }
        if (pMin != null) {
            min = pMin;
            doMin = true;
        }

        HashMap<String, Double> regionMap = CoverageUtilities.getRegionParamsFromGridCoverage(inRaster);
        int nCols = regionMap.get(CoverageUtilities.COLS).intValue();
        int nRows = regionMap.get(CoverageUtilities.ROWS).intValue();

        RenderedImage geodataRI = inRaster.getRenderedImage();
        RandomIter geodataIter = RandomIterFactory.create(geodataRI, null);

        if (inMask != null) {
            RenderedImage maskRI = inMask.getRenderedImage();
            maskIter = RandomIterFactory.create(maskRI, null);
        }

        WritableRaster outWR = CoverageUtilities.renderedImage2WritableRaster(geodataRI, false);
        WritableRandomIter outIter = RandomIterFactory.createWritable(outWR, null);

        pm.beginTask("Processing map...", nRows);
        for( int i = 0; i < nRows; i++ ) {
            if (isCanceled(pm)) {
                return;
            }
            for( int j = 0; j < nCols; j++ ) {
                if (maskIter != null) {
                    double maskValue = maskIter.getSampleDouble(j, i, 0);
                    if (!doInverse) {
                        if (isNovalue(maskValue)) {
                            outIter.setSample(j, i, 0, JGTConstants.doubleNovalue);
                            continue;
                        }
                    } else {
                        if (!isNovalue(maskValue)) {
                            outIter.setSample(j, i, 0, JGTConstants.doubleNovalue);
                            continue;
                        }
                    }
                }
                double value = geodataIter.getSampleDouble(j, i, 0);
                if (doMax && value > max) {
                    outIter.setSample(j, i, 0, JGTConstants.doubleNovalue);
                    continue;
                }
                if (doMin && value < min) {
                    outIter.setSample(j, i, 0, JGTConstants.doubleNovalue);
                    continue;
                }
            }
            pm.worked(1);
        }
        pm.done();

        outRaster = CoverageUtilities.buildCoverage("pitfiller", outWR, regionMap, inRaster.getCoordinateReferenceSystem());

    }
}
