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

import static org.jgrasstools.gears.libs.modules.Variables.TCA;
import static org.jgrasstools.gears.libs.modules.Variables.TCA_CONVERGENT;
import static org.jgrasstools.gears.libs.modules.Variables.TCA_SLOPE;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_AUTHORCONTACTS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_AUTHORNAMES;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_KEYWORDS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_LABEL;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_LICENSE;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_NAME;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_STATUS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_inFlow_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_inSlope_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_inTc3_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_inTca_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_outNet_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_pExp_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_pMode_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSEXTRACTNETWORK_pThres_DESCRIPTION;
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
import org.jgrasstools.hortonmachine.modules.network.extractnetwork.OmsExtractNetwork;

@Description(OMSEXTRACTNETWORK_DESCRIPTION)
@Author(name = OMSEXTRACTNETWORK_AUTHORNAMES, contact = OMSEXTRACTNETWORK_AUTHORCONTACTS)
@Keywords(OMSEXTRACTNETWORK_KEYWORDS)
@Label(OMSEXTRACTNETWORK_LABEL)
@Name("_" + OMSEXTRACTNETWORK_NAME)
@Status(OMSEXTRACTNETWORK_STATUS)
@License(OMSEXTRACTNETWORK_LICENSE)
public class ExtractNetwork extends JGTModel {

    @Description(OMSEXTRACTNETWORK_inTca_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inTca = null;

    @Description(OMSEXTRACTNETWORK_inFlow_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inFlow = null;

    @Description(OMSEXTRACTNETWORK_inSlope_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inSlope = null;

    @Description(OMSEXTRACTNETWORK_inTc3_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inTc3 = null;

    @Description(OMSEXTRACTNETWORK_pThres_DESCRIPTION)
    @In
    public double pThres = 0;

    @Description(OMSEXTRACTNETWORK_pMode_DESCRIPTION)
    @UI("combo:" + TCA + "," + TCA_SLOPE + "," + TCA_CONVERGENT)
    @In
    public String pMode = TCA;

    @Description(OMSEXTRACTNETWORK_pExp_DESCRIPTION)
    @In
    public double pExp = 0.5;

    @Description(OMSEXTRACTNETWORK_outNet_DESCRIPTION)
    @UI(JGTConstants.FILEOUT_UI_HINT)
    @In
    public String outNet = null;

    @Execute
    public void process() throws Exception {
        OmsExtractNetwork extractnetwork = new OmsExtractNetwork();
        extractnetwork.inTca = getRaster(inTca);
        extractnetwork.inFlow = getRaster(inFlow);
        extractnetwork.inSlope = getRaster(inSlope);
        extractnetwork.inTc3 = getRaster(inTc3);
        extractnetwork.pThres = pThres;
        extractnetwork.pMode = pMode;
        extractnetwork.pExp = pExp;
        extractnetwork.pm = pm;
        extractnetwork.doProcess = doProcess;
        extractnetwork.doReset = doReset;
        extractnetwork.process();
        dumpRaster(extractnetwork.outNet, outNet);
    }
    

}
