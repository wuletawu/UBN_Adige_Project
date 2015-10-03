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

import static org.jgrasstools.gears.i18n.GearsMessages.OMSINTERSECTIONFINDER_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSINTERSECTIONFINDER_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSINTERSECTIONFINDER_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSINTERSECTIONFINDER_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSINTERSECTIONFINDER_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSINTERSECTIONFINDER_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSINTERSECTIONFINDER_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSINTERSECTIONFINDER_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSINTERSECTIONFINDER_inMap_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSINTERSECTIONFINDER_outLinesMap_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSINTERSECTIONFINDER_outPointsMap_DESCRIPTION;
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
import org.jgrasstools.gears.modules.v.intersections.OmsIntersectionFinder;

@Description(OMSINTERSECTIONFINDER_DESCRIPTION)
@Author(name = OMSINTERSECTIONFINDER_AUTHORNAMES, contact = OMSINTERSECTIONFINDER_AUTHORCONTACTS)
@Keywords(OMSINTERSECTIONFINDER_KEYWORDS)
@Label(OMSINTERSECTIONFINDER_LABEL)
@Name("_" + OMSINTERSECTIONFINDER_NAME)
@Status(OMSINTERSECTIONFINDER_STATUS)
@License(OMSINTERSECTIONFINDER_LICENSE)
public class IntersectionFinder extends JGTModel {

    @Description(OMSINTERSECTIONFINDER_inMap_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inMap = null;

    @Description(OMSINTERSECTIONFINDER_outPointsMap_DESCRIPTION)
    @UI(JGTConstants.FILEOUT_UI_HINT)
    @In
    public String outPointsMap = null;

    @Description(OMSINTERSECTIONFINDER_outLinesMap_DESCRIPTION)
    @UI(JGTConstants.FILEOUT_UI_HINT)
    @In
    public String outLinesMap = null;

    @Execute
    public void process() throws Exception {
        OmsIntersectionFinder intersectionfinder = new OmsIntersectionFinder();
        intersectionfinder.inMap = getVector(inMap);
        intersectionfinder.pm = pm;
        intersectionfinder.doProcess = doProcess;
        intersectionfinder.doReset = doReset;
        intersectionfinder.process();
        dumpVector(intersectionfinder.outPointsMap, outPointsMap);
        dumpVector(intersectionfinder.outLinesMap, outLinesMap);
    }
}
