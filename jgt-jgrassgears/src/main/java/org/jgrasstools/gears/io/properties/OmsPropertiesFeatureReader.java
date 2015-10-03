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
package org.jgrasstools.gears.io.properties;

import static org.jgrasstools.gears.i18n.GearsMessages.OMSPROPERTIESFEATUREREADER_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPROPERTIESFEATUREREADER_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPROPERTIESFEATUREREADER_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPROPERTIESFEATUREREADER_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPROPERTIESFEATUREREADER_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPROPERTIESFEATUREREADER_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPROPERTIESFEATUREREADER_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPROPERTIESFEATUREREADER_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPROPERTIESFEATUREREADER_UI;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPROPERTIESFEATUREREADER_file_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSPROPERTIESFEATUREREADER_geodata_DESCRIPTION;

import java.io.File;
import java.io.IOException;

import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;
import oms3.annotations.UI;

import org.geotools.data.FeatureReader;
import org.geotools.data.Query;
import org.geotools.data.Transaction;
import org.geotools.data.property.PropertyDataStore;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.DefaultFeatureCollection;
import org.geotools.feature.FeatureCollection;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.utils.files.FileUtilities;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

@Description(OMSPROPERTIESFEATUREREADER_DESCRIPTION)
@Author(name = OMSPROPERTIESFEATUREREADER_AUTHORNAMES, contact = OMSPROPERTIESFEATUREREADER_AUTHORCONTACTS)
@Keywords(OMSPROPERTIESFEATUREREADER_KEYWORDS)
@Label(OMSPROPERTIESFEATUREREADER_LABEL)
@Name(OMSPROPERTIESFEATUREREADER_NAME)
@Status(OMSPROPERTIESFEATUREREADER_STATUS)
@License(OMSPROPERTIESFEATUREREADER_LICENSE)
@UI(OMSPROPERTIESFEATUREREADER_UI)
public class OmsPropertiesFeatureReader extends JGTModel {

    @Description(OMSPROPERTIESFEATUREREADER_file_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String file = null;

    @Description(OMSPROPERTIESFEATUREREADER_geodata_DESCRIPTION)
    @Out
    public SimpleFeatureCollection geodata = null;

    @Execute
    public void readFeatureCollection() throws IOException {
        if (!concatOr(geodata == null, doReset)) {
            return;
        }

        /*
         * Works on types:
         * 
         * _=id:Integer,name:String,geom:Point
         * fid1=1|jody garnett|POINT(0 0)
         * fid2=2|brent|POINT(10 10)
         * fid3=3|dave|POINT(20 20)
         * fid4=4|justin deolivera|POINT(30 30)
         */

        FeatureReader<SimpleFeatureType, SimpleFeature> reader = null;
        try {
            File propertiesFile = new File(file);
            pm.beginTask("Reading features from properties file: " + propertiesFile.getName(), -1);
            PropertyDataStore store = new PropertyDataStore(propertiesFile.getParentFile());

            String name = FileUtilities.getNameWithoutExtention(propertiesFile);

            geodata = new DefaultFeatureCollection();
            Query query = new Query(name);
            reader = store.getFeatureReader(query, Transaction.AUTO_COMMIT);
            while( reader.hasNext() ) {
                SimpleFeature feature = reader.next();
                ((DefaultFeatureCollection) geodata).add(feature);
            }
        } finally {
            pm.done();
            if (reader != null)
                reader.close();
        }
    }

    /**
     * Fast read access mode. 
     * 
     * @param path the properties file path.
     * @return the read {@link FeatureCollection}.
     * @throws IOException
     */
    public static SimpleFeatureCollection readPropertiesfile( String path ) throws IOException {

        OmsPropertiesFeatureReader reader = new OmsPropertiesFeatureReader();
        reader.file = path;
        reader.readFeatureCollection();

        return reader.geodata;
    }

}
