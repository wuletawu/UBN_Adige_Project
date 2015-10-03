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
package org.jgrasstools.hortonmachine.modules.geomorphology.aspect;

import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSASPECT_AUTHORCONTACTS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSASPECT_AUTHORNAMES;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSASPECT_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSASPECT_DOCUMENTATION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSASPECT_KEYWORDS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSASPECT_LABEL;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSASPECT_LICENSE;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSASPECT_NAME;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSASPECT_STATUS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSASPECT_doRadiants_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSASPECT_doRound_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSASPECT_inElev_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSASPECT_outAspect_DESCRIPTION;

import java.io.File;

import javax.media.jai.iterator.RandomIter;
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

import org.jgrasstools.gears.libs.modules.GridNode;
import org.jgrasstools.gears.libs.modules.JGTModelIM;
import org.jgrasstools.gears.utils.colors.ColorTables;
import org.jgrasstools.gears.utils.math.NumericsUtilities;

@Description(OMSASPECT_DESCRIPTION)
@Documentation(OMSASPECT_DOCUMENTATION)
@Author(name = OMSASPECT_AUTHORNAMES, contact = OMSASPECT_AUTHORCONTACTS)
@Keywords(OMSASPECT_KEYWORDS)
@Label(OMSASPECT_LABEL)
@Name(OMSASPECT_NAME)
@Status(OMSASPECT_STATUS)
@License(OMSASPECT_LICENSE)
public class OmsAspectIM extends JGTModelIM {
    @Description(OMSASPECT_inElev_DESCRIPTION)
    @In
    public String inElev = null;

    @Description(OMSASPECT_doRadiants_DESCRIPTION)
    @In
    public boolean doRadiants = false;

    @Description(OMSASPECT_doRound_DESCRIPTION)
    @In
    public boolean doRound = false;

    @Description(OMSASPECT_outAspect_DESCRIPTION)
    @Out
    public String outAspect = null;

    private double radtodeg;

    @Execute
    public void process() throws Exception {
        checkNull(inElev);
        radtodeg = NumericsUtilities.RADTODEG;
        if (doRadiants) {
            radtodeg = 1.0;
        }

        cellBuffer = 1;
        addSource(new File(inElev));
        addDestination(new File(outAspect));

        processByTileCells();

        makeMosaic();
        makeStyle(ColorTables.aspect, 0, 360);
        
        dispose();
    }

    @Override
    protected void processCell( int readCol, int readRow, int writeCol, int writeRow, int readCols, int readRows, int writeCols,
            int writeRows ) {

        RandomIter elevIter = inRasterIterators.get(0);
        GridNode node = new GridNode(elevIter, readCols, readRows, xRes, yRes, readCol, readRow);
        double aspect = OmsAspect.calculateAspect(node, radtodeg, doRound);
        WritableRandomIter outDataIter = outRasters.get(0);
        outDataIter.setSample(writeCol, writeRow, 0, aspect);
    }

}
