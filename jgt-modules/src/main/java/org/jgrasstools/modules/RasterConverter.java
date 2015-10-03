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

import static org.jgrasstools.gears.i18n.GearsMessages.GENERIC_pCols_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.GENERIC_pEast_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.GENERIC_pNorth_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.GENERIC_pRows_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.GENERIC_pSouth_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.GENERIC_pWest_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.GENERIC_pXres_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.GENERIC_pYres_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERCONVERTER_inRaster_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERCONVERTER_outRaster_DESCRIPTION;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Name;
import oms3.annotations.UI;

import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.modules.r.rasterconverter.OmsRasterConverter;

@Name("rconvert")
public class RasterConverter extends OmsRasterConverter {

    @Description(OMSRASTERCONVERTER_inRaster_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inRaster;

    @Description(GENERIC_pNorth_DESCRIPTION)
    @UI(JGTConstants.PROCESS_NORTH_UI_HINT)
    @In
    public Double pNorth = null;

    @Description(GENERIC_pSouth_DESCRIPTION)
    @UI(JGTConstants.PROCESS_SOUTH_UI_HINT)
    @In
    public Double pSouth = null;

    @Description(GENERIC_pWest_DESCRIPTION)
    @UI(JGTConstants.PROCESS_WEST_UI_HINT)
    @In
    public Double pWest = null;

    @Description(GENERIC_pEast_DESCRIPTION)
    @UI(JGTConstants.PROCESS_EAST_UI_HINT)
    @In
    public Double pEast = null;
    
    @Description(GENERIC_pXres_DESCRIPTION)
    @UI(JGTConstants.PROCESS_XRES_UI_HINT)
    @In
    public Double pXres = null;

    @Description(GENERIC_pYres_DESCRIPTION)
    @UI(JGTConstants.PROCESS_YRES_UI_HINT)
    @In
    public Double pYres = null;

    @Description(GENERIC_pRows_DESCRIPTION)
    @UI(JGTConstants.PROCESS_ROWS_UI_HINT)
    @In
    public Integer pRows = null;

    @Description(GENERIC_pCols_DESCRIPTION)
    @UI(JGTConstants.PROCESS_COLS_UI_HINT)
    @In
    public Integer pCols = null;

    @Description(OMSRASTERCONVERTER_outRaster_DESCRIPTION)
    @UI(JGTConstants.FILEOUT_UI_HINT)
    @In
    public String outRaster;

    @Execute
    public void process() throws Exception {
        
        OmsRasterReader rasterreader = new OmsRasterReader();
        rasterreader.file = inRaster;
        rasterreader.pNorth = pNorth;
        rasterreader.pSouth = pSouth;
        rasterreader.pWest = pWest;
        rasterreader.pEast = pEast;
        rasterreader.pXres = pXres;
        rasterreader.pYres = pYres;
        rasterreader.pRows = pRows;
        rasterreader.pCols = pCols;
        rasterreader.pm = pm;
        rasterreader.process();
        OmsRasterConverter rasterconverter = new OmsRasterConverter();
        rasterconverter.inRaster = rasterreader.outRaster;
        rasterconverter.pm = pm;
        rasterconverter.doProcess = doProcess;
        rasterconverter.doReset = doReset;
        rasterconverter.process();
        dumpRaster(rasterconverter.outRaster, outRaster);
    }
}
