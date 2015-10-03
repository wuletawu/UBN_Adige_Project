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
package org.jgrasstools.gears.io.vectorwriter;

import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORWRITER_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORWRITER_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORWRITER_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORWRITER_DOCUMENTATION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORWRITER_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORWRITER_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORWRITER_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORWRITER_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORWRITER_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORWRITER_file_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORWRITER_inVector_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORWRITER_pType_DESCRIPTION;

import java.io.File;
import java.io.IOException;

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

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureCollection;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureWriter;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;

@Description(OMSVECTORWRITER_DESCRIPTION)
@Documentation(OMSVECTORWRITER_DOCUMENTATION)
@Author(name = OMSVECTORWRITER_AUTHORNAMES, contact = OMSVECTORWRITER_AUTHORCONTACTS)
@Keywords(OMSVECTORWRITER_KEYWORDS)
@Label(OMSVECTORWRITER_LABEL)
@Name(OMSVECTORWRITER_NAME)
@Status(OMSVECTORWRITER_STATUS)
@License(OMSVECTORWRITER_LICENSE)
public class OmsVectorWriter extends JGTModel {

    @Description(OMSVECTORWRITER_inVector_DESCRIPTION)
    @In
    public SimpleFeatureCollection inVector = null;

    @Description(OMSVECTORWRITER_pType_DESCRIPTION)
    @In
    // currently not used, for future compatibility
    public String pType = null;

    @Description(OMSVECTORWRITER_file_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String file = null;

    @Execute
    public void process() throws IOException {
        checkNull(file);

        File vectorFile = new File(file);
        if (inVector.size() == 0) {
            pm.message("Warning, not writing an empty vector to file: " + vectorFile.getName());
            return;
        }
        String name = vectorFile.getName();
        if (name.toLowerCase().endsWith("shp") || (pType != null && pType.equals(JGTConstants.SHP))) {
            OmsShapefileFeatureWriter.writeShapefile(vectorFile.getAbsolutePath(), inVector, pm);
        } else {
            throw new IOException("Format is currently not supported for file: " + name);
        }
    }

    /**
     * Fast write access mode. 
     * 
     * @param path the vector file path.
     * @param featureCollection the {@link FeatureCollection} to write.
     * @throws IOException
     */
    public static void writeVector( String path, SimpleFeatureCollection featureCollection ) throws IOException {
        OmsVectorWriter writer = new OmsVectorWriter();
        writer.file = path;
        writer.inVector = featureCollection;
        writer.process();
    }

}
