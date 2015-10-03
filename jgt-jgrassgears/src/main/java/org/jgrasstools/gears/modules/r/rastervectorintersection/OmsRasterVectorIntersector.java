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
package org.jgrasstools.gears.modules.r.rastervectorintersection;

import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERVECTORINTERSECTOR_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERVECTORINTERSECTOR_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERVECTORINTERSECTOR_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERVECTORINTERSECTOR_DOCUMENTATION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERVECTORINTERSECTOR_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERVECTORINTERSECTOR_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERVECTORINTERSECTOR_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERVECTORINTERSECTOR_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERVECTORINTERSECTOR_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERVECTORINTERSECTOR_doInverse_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERVECTORINTERSECTOR_inRaster_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERVECTORINTERSECTOR_inVector_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERVECTORINTERSECTOR_outRaster_DESCRIPTION;
import static org.jgrasstools.gears.utils.geometry.GeometryUtilities.getGeometryType;
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
import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.libs.exceptions.ModelsRuntimeException;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.modules.r.cutout.OmsCutOut;
import org.jgrasstools.gears.modules.r.scanline.OmsScanLineRasterizer;
import org.jgrasstools.gears.utils.RegionMap;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.feature.type.GeometryType;

@Description(OMSRASTERVECTORINTERSECTOR_DESCRIPTION)
@Documentation(OMSRASTERVECTORINTERSECTOR_DOCUMENTATION)
@Author(name = OMSRASTERVECTORINTERSECTOR_AUTHORNAMES, contact = OMSRASTERVECTORINTERSECTOR_AUTHORCONTACTS)
@Keywords(OMSRASTERVECTORINTERSECTOR_KEYWORDS)
@Label(OMSRASTERVECTORINTERSECTOR_LABEL)
@Name(OMSRASTERVECTORINTERSECTOR_NAME)
@Status(OMSRASTERVECTORINTERSECTOR_STATUS)
@License(OMSRASTERVECTORINTERSECTOR_LICENSE)
public class OmsRasterVectorIntersector extends JGTModel {

    @Description(OMSRASTERVECTORINTERSECTOR_inVector_DESCRIPTION)
    @In
    public SimpleFeatureCollection inVector = null;

    @Description(OMSRASTERVECTORINTERSECTOR_inRaster_DESCRIPTION)
    @In
    public GridCoverage2D inRaster;

    @Description(OMSRASTERVECTORINTERSECTOR_doInverse_DESCRIPTION)
    @In
    public boolean doInverse = false;

    @Description(OMSRASTERVECTORINTERSECTOR_outRaster_DESCRIPTION)
    @Out
    public GridCoverage2D outRaster;

    @Execute
    public void process() throws Exception {
        checkNull(inRaster, inVector);

        SimpleFeatureType schema = inVector.getSchema();
        GeometryType type = schema.getGeometryDescriptor().getType();
        if (getGeometryType(type) != org.jgrasstools.gears.utils.geometry.GeometryType.POLYGON && getGeometryType(type) != org.jgrasstools.gears.utils.geometry.GeometryType.MULTIPOLYGON) {
            throw new ModelsRuntimeException("The module works only with polygon vectors.", this);
        }

        RegionMap regionMap = CoverageUtilities.getRegionParamsFromGridCoverage(inRaster);
        OmsScanLineRasterizer raster = new OmsScanLineRasterizer();
        raster.inVector = inVector;
        raster.pCols = regionMap.getCols();
        raster.pRows = regionMap.getRows();
        raster.pNorth = regionMap.getNorth();
        raster.pSouth = regionMap.getSouth();
        raster.pEast = regionMap.getEast();
        raster.pWest = regionMap.getWest();
        raster.pValue = 1.0;
        raster.process();
        GridCoverage2D rasterizedVector = raster.outRaster;

        OmsCutOut cutout = new OmsCutOut();
        cutout.pm = pm;
        cutout.inRaster = inRaster;
        cutout.inMask = rasterizedVector;
        cutout.doInverse = doInverse;
        cutout.process();
        outRaster = cutout.outRaster;
    }
}
