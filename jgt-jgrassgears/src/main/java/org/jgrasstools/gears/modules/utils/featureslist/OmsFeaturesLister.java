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
package org.jgrasstools.gears.modules.utils.featureslist;

import static org.jgrasstools.gears.i18n.GearsMessages.OMSFEATURESLISTER_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSFEATURESLISTER_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSFEATURESLISTER_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSFEATURESLISTER_DOCUMENTATION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSFEATURESLISTER_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSFEATURESLISTER_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSFEATURESLISTER_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSFEATURESLISTER_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSFEATURESLISTER_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSFEATURESLISTER_inFiles_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSFEATURESLISTER_outFC_DESCRIPTION;

import java.util.ArrayList;
import java.util.List;

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
import oms3.annotations.UI;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.vectorreader.OmsVectorReader;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;

@Description(OMSFEATURESLISTER_DESCRIPTION)
@Documentation(OMSFEATURESLISTER_DOCUMENTATION)
@Author(name = OMSFEATURESLISTER_AUTHORNAMES, contact = OMSFEATURESLISTER_AUTHORCONTACTS)
@Keywords(OMSFEATURESLISTER_KEYWORDS)
@Label(OMSFEATURESLISTER_LABEL)
@Name(OMSFEATURESLISTER_NAME)
@Status(OMSFEATURESLISTER_STATUS)
@License(OMSFEATURESLISTER_LICENSE)
public class OmsFeaturesLister extends JGTModel {

    @Description(OMSFEATURESLISTER_inFiles_DESCRIPTION)
    @UI(JGTConstants.FILESPATHLIST_UI_HINT)
    @In
    public List<String> inFiles;

    @Description(OMSFEATURESLISTER_outFC_DESCRIPTION)
    @Out
    public List<SimpleFeatureCollection> outFC = null;

    @Execute
    public void process() throws Exception {

        outFC = new ArrayList<SimpleFeatureCollection>();

        for( String file : inFiles ) {
            SimpleFeatureCollection featureCollection = OmsVectorReader.readVector(file);
            outFC.add(featureCollection);
        }

    }

}
