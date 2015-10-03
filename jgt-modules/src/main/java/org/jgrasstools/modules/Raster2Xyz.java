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

import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTER2XYZ_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTER2XYZ_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTER2XYZ_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTER2XYZ_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTER2XYZ_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTER2XYZ_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTER2XYZ_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTER2XYZ_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTER2XYZ_doRemovenv_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTER2XYZ_inFile_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTER2XYZ_inRaster_DESCRIPTION;
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

import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.modules.r.raster2xyz.OmsRaster2Xyz;

@Description(OMSRASTER2XYZ_DESCRIPTION)
@Author(name = OMSRASTER2XYZ_AUTHORNAMES, contact = OMSRASTER2XYZ_AUTHORCONTACTS)
@Keywords(OMSRASTER2XYZ_KEYWORDS)
@Label(OMSRASTER2XYZ_LABEL)
@Name("_" + OMSRASTER2XYZ_NAME)
@Status(OMSRASTER2XYZ_STATUS)
@License(OMSRASTER2XYZ_LICENSE)
public class Raster2Xyz extends JGTModel {

    @Description(OMSRASTER2XYZ_inRaster_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inRaster;

    @Description(OMSRASTER2XYZ_inFile_DESCRIPTION)
    @UI(JGTConstants.FILEOUT_UI_HINT)
    @In
    public String inFile;

    @Description(OMSRASTER2XYZ_doRemovenv_DESCRIPTION)
    @In
    public boolean doRemovenv = true;

    @SuppressWarnings("nls")
    @Execute
    public void process() throws Exception {
        OmsRaster2Xyz raster2xyz = new OmsRaster2Xyz();
        raster2xyz.inRaster = getRaster(inRaster);
        raster2xyz.inFile = inFile;
        raster2xyz.doRemovenv = doRemovenv;
        raster2xyz.pm = pm;
        raster2xyz.doProcess = doProcess;
        raster2xyz.doReset = doReset;
        raster2xyz.process();
    }
}
