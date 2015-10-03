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

import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERNULL_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERNULL_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERNULL_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERNULL_DOCUMENTATION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERNULL_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERNULL_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERNULL_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERNULL_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERNULL_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERNULL_inRaster_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERNULL_outRaster_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERNULL_pNull_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERNULL_pValue_DESCRIPTION;
import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
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
import org.jgrasstools.gears.modules.r.rasternull.OmsRasterNull;

@Description(OMSRASTERNULL_DESCRIPTION)
@Documentation(OMSRASTERNULL_DOCUMENTATION)
@Author(name = OMSRASTERNULL_AUTHORNAMES, contact = OMSRASTERNULL_AUTHORCONTACTS)
@Keywords(OMSRASTERNULL_KEYWORDS)
@Label(OMSRASTERNULL_LABEL)
@Name(OMSRASTERNULL_NAME)
@Status(OMSRASTERNULL_STATUS)
@License(OMSRASTERNULL_LICENSE)
public class RasterNull extends JGTModel {

    @Description(OMSRASTERNULL_inRaster_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inRaster;

    @Description(OMSRASTERNULL_pValue_DESCRIPTION)
    @In
    public Double pValue = null;

    @Description(OMSRASTERNULL_pNull_DESCRIPTION)
    @In
    public Double pNull = null;

    @Description(OMSRASTERNULL_outRaster_DESCRIPTION)
    @UI(JGTConstants.FILEOUT_UI_HINT)
    @In
    public String outRaster;

    @Execute
    public void process() throws Exception {
        OmsRasterNull rasternull = new OmsRasterNull();
        rasternull.inRaster = getRaster(inRaster);
        rasternull.pValue = pValue;
        rasternull.pNull = pNull;
        rasternull.pm = pm;
        rasternull.doProcess = doProcess;
        rasternull.doReset = doReset;
        rasternull.process();
        dumpRaster(rasternull.outRaster, outRaster);
    }

}
