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

import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_inMap10_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_inMap11_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_inMap12_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_inMap1_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_inMap2_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_inMap3_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_inMap4_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_inMap5_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_inMap6_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_inMap7_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_inMap8_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_inMap9_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_outMap_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMOSAIC12_pInterpolation_DESCRIPTION;
import static org.jgrasstools.gears.libs.modules.Variables.BICUBIC;
import static org.jgrasstools.gears.libs.modules.Variables.BILINEAR;
import static org.jgrasstools.gears.libs.modules.Variables.NEAREST_NEIGHTBOUR;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

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

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.io.rasterwriter.OmsRasterWriter;
import org.jgrasstools.gears.libs.exceptions.ModelsIllegalargumentException;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.modules.r.mosaic.OmsMosaic;

@Description(OMSMOSAIC12_DESCRIPTION)
@Author(name = OMSMOSAIC12_AUTHORNAMES, contact = OMSMOSAIC12_AUTHORCONTACTS)
@Keywords(OMSMOSAIC12_KEYWORDS)
@Label(OMSMOSAIC12_LABEL)
@Name("_" + OMSMOSAIC12_NAME)
@Status(OMSMOSAIC12_STATUS)
@License(OMSMOSAIC12_LICENSE)
public class Mosaic12 extends JGTModel {

    @Description(OMSMOSAIC12_inMap1_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inMap1;

    @Description(OMSMOSAIC12_inMap2_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inMap2;

    @Description(OMSMOSAIC12_inMap3_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inMap3;

    @Description(OMSMOSAIC12_inMap4_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inMap4;

    @Description(OMSMOSAIC12_inMap5_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inMap5;

    @Description(OMSMOSAIC12_inMap6_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inMap6;

    @Description(OMSMOSAIC12_inMap7_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inMap7;

    @Description(OMSMOSAIC12_inMap8_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inMap8;

    @Description(OMSMOSAIC12_inMap9_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inMap9;

    @Description(OMSMOSAIC12_inMap10_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inMap10;

    @Description(OMSMOSAIC12_inMap11_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inMap11;

    @Description(OMSMOSAIC12_inMap12_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String inMap12;

    @Description(OMSMOSAIC12_pInterpolation_DESCRIPTION)
    @UI("combo:" + NEAREST_NEIGHTBOUR + "," + BILINEAR + "," + BICUBIC)
    @In
    public String pInterpolation = NEAREST_NEIGHTBOUR;

    @Description(OMSMOSAIC12_outMap_DESCRIPTION)
    @UI(JGTConstants.FILEOUT_UI_HINT)
    @In
    public String outMap = null;

    public GridCoverage2D outRaster;

    @Execute
    public void process() throws Exception {
        checkNull(outMap);

        List<File> filesList = new ArrayList<File>();
        checkMap(filesList, inMap1);
        checkMap(filesList, inMap2);
        checkMap(filesList, inMap3);
        checkMap(filesList, inMap4);
        checkMap(filesList, inMap5);
        checkMap(filesList, inMap6);
        checkMap(filesList, inMap7);
        checkMap(filesList, inMap8);
        checkMap(filesList, inMap9);
        checkMap(filesList, inMap10);
        checkMap(filesList, inMap11);
        checkMap(filesList, inMap12);

        if (filesList.size() < 2) {
            throw new ModelsIllegalargumentException("The patching module needs at least two maps to be patched.", this, pm);
        }

        OmsMosaic mosaic = new OmsMosaic();
        mosaic.inFiles = filesList;
        mosaic.pm = pm;
        mosaic.process();

        outRaster = mosaic.outRaster;
        OmsRasterWriter.writeRaster(outMap, outRaster);
    }

    private void checkMap( List<File> filesList, String inMap ) {
        if (inMap != null) {
            File tmp = new File(inMap);
            if (tmp.exists()) {
                filesList.add(tmp);
            }
        }
    }

}
