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
package org.jgrasstools.modules;

import static org.jgrasstools.gears.i18n.GearsMessages.OMSGRIDGEOMETRYREADER_pCode_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSGRIDGEOMETRYREADER_pEast_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSGRIDGEOMETRYREADER_pNorth_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSGRIDGEOMETRYREADER_pSouth_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSGRIDGEOMETRYREADER_pWest_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSGRIDGEOMETRYREADER_pXres_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSGRIDGEOMETRYREADER_pYres_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPOINTSRASTERIZER_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPOINTSRASTERIZER_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPOINTSRASTERIZER_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPOINTSRASTERIZER_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPOINTSRASTERIZER_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPOINTSRASTERIZER_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPOINTSRASTERIZER_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPOINTSRASTERIZER_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPOINTSRASTERIZER_fCat_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPOINTSRASTERIZER_inVector_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPOINTSRASTERIZER_outRaster_DESCRIPTION;
import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Status;
import oms3.annotations.UI;

import org.geotools.coverage.grid.GridGeometry2D;
import org.jgrasstools.gears.io.gridgeometryreader.OmsGridGeometryReader;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.modules.r.pointsrasterizer.OmsPointsRasterizer;

@Description(OMSPOINTSRASTERIZER_DESCRIPTION)
@Author(name = OMSPOINTSRASTERIZER_AUTHORNAMES, contact = OMSPOINTSRASTERIZER_AUTHORCONTACTS)
@Keywords(OMSPOINTSRASTERIZER_KEYWORDS)
@Label(OMSPOINTSRASTERIZER_LABEL)
@Name("_" + OMSPOINTSRASTERIZER_NAME)
@Status(OMSPOINTSRASTERIZER_STATUS)
@License(OMSPOINTSRASTERIZER_LICENSE)
public class PointsRasterizer extends JGTModel {

    @Description(OMSPOINTSRASTERIZER_inVector_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inVector = null;

    @Description(OMSGRIDGEOMETRYREADER_pNorth_DESCRIPTION)
    @UI(JGTConstants.PROCESS_NORTH_UI_HINT)
    @In
    public Double pNorth = null;

    @Description(OMSGRIDGEOMETRYREADER_pSouth_DESCRIPTION)
    @UI(JGTConstants.PROCESS_SOUTH_UI_HINT)
    @In
    public Double pSouth = null;

    @Description(OMSGRIDGEOMETRYREADER_pWest_DESCRIPTION)
    @UI(JGTConstants.PROCESS_WEST_UI_HINT)
    @In
    public Double pWest = null;

    @Description(OMSGRIDGEOMETRYREADER_pEast_DESCRIPTION)
    @UI(JGTConstants.PROCESS_EAST_UI_HINT)
    @In
    public Double pEast = null;

    @Description(OMSGRIDGEOMETRYREADER_pXres_DESCRIPTION)
    @UI(JGTConstants.PROCESS_XRES_UI_HINT)
    @In
    public Double pXres = null;

    @Description(OMSGRIDGEOMETRYREADER_pYres_DESCRIPTION)
    @UI(JGTConstants.PROCESS_YRES_UI_HINT)
    @In
    public Double pYres = null;

    @Description(OMSGRIDGEOMETRYREADER_pCode_DESCRIPTION)
    @UI(JGTConstants.CRS_UI_HINT)
    @In
    public String pCode;

    @Description(OMSPOINTSRASTERIZER_fCat_DESCRIPTION)
    @In
    public String fCat;

    @Description(OMSPOINTSRASTERIZER_outRaster_DESCRIPTION)
    @UI(JGTConstants.FILEOUT_UI_HINT)
    @In
    public String outRaster;

    @Execute
    public void process() throws Exception {
        OmsGridGeometryReader gridgeometryreader = new OmsGridGeometryReader();
        gridgeometryreader.pNorth = pNorth;
        gridgeometryreader.pSouth = pSouth;
        gridgeometryreader.pWest = pWest;
        gridgeometryreader.pEast = pEast;
        gridgeometryreader.pXres = pXres;
        gridgeometryreader.pYres = pYres;
        gridgeometryreader.pCode = pCode;
        gridgeometryreader.pm = pm;
        gridgeometryreader.doProcess = doProcess;
        gridgeometryreader.doReset = doReset;
        gridgeometryreader.process();
        GridGeometry2D outGridgeom = gridgeometryreader.outGridgeom;

        OmsPointsRasterizer pointsrasterizer = new OmsPointsRasterizer();
        pointsrasterizer.inVector = getVector(inVector);
        pointsrasterizer.inGrid = outGridgeom;
        pointsrasterizer.fCat = fCat;
        pointsrasterizer.pm = pm;
        pointsrasterizer.doProcess = doProcess;
        pointsrasterizer.doReset = doReset;
        pointsrasterizer.process();
        dumpRaster(pointsrasterizer.outRaster, outRaster);
    }
}
