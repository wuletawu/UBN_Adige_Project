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

import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSMARKOUTLETS_AUTHORCONTACTS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSMARKOUTLETS_AUTHORNAMES;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSMARKOUTLETS_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSMARKOUTLETS_KEYWORDS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSMARKOUTLETS_LABEL;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSMARKOUTLETS_LICENSE;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSMARKOUTLETS_NAME;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSMARKOUTLETS_STATUS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSMARKOUTLETS_inFlow_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSMARKOUTLETS_outFlow_DESCRIPTION;
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
import org.jgrasstools.hortonmachine.modules.demmanipulation.markoutlets.OmsMarkoutlets;

@Description(OMSMARKOUTLETS_DESCRIPTION)
@Author(name = OMSMARKOUTLETS_AUTHORNAMES, contact = OMSMARKOUTLETS_AUTHORCONTACTS)
@Keywords(OMSMARKOUTLETS_KEYWORDS)
@Label(OMSMARKOUTLETS_LABEL)
@Name("_" + OMSMARKOUTLETS_NAME)
@Status(OMSMARKOUTLETS_STATUS)
@License(OMSMARKOUTLETS_LICENSE)
public class Markoutlets extends JGTModel {
    @Description(OMSMARKOUTLETS_inFlow_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inFlow = null;

    @Description(OMSMARKOUTLETS_outFlow_DESCRIPTION)
    @UI(JGTConstants.FILEOUT_UI_HINT)
    @In
    public String outFlow = null;

    @Execute
    public void process() throws Exception {
        OmsMarkoutlets omsmarkoutlets = new OmsMarkoutlets();
        omsmarkoutlets.inFlow = getRaster(inFlow);
        omsmarkoutlets.pm = pm;
        omsmarkoutlets.doProcess = doProcess;
        omsmarkoutlets.doReset = doReset;
        omsmarkoutlets.process();
        dumpRaster(omsmarkoutlets.outFlow, outFlow);
    }
}
