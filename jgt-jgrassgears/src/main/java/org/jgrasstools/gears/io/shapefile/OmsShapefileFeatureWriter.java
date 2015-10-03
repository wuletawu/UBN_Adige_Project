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
package org.jgrasstools.gears.io.shapefile;

import static org.jgrasstools.gears.i18n.GearsMessages.OMSSHAPEFILEFEATUREWRITER_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSSHAPEFILEFEATUREWRITER_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSSHAPEFILEFEATUREWRITER_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSSHAPEFILEFEATUREWRITER_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSSHAPEFILEFEATUREWRITER_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSSHAPEFILEFEATUREWRITER_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSSHAPEFILEFEATUREWRITER_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSSHAPEFILEFEATUREWRITER_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSSHAPEFILEFEATUREWRITER_UI;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSSHAPEFILEFEATUREWRITER_doIndex_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSSHAPEFILEFEATUREWRITER_file_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSSHAPEFILEFEATUREWRITER_geodata_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSSHAPEFILEFEATUREWRITER_pType_DESCRIPTION;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Map;

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

import org.geotools.data.DataStore;
import org.geotools.data.DefaultTransaction;
import org.geotools.data.FileDataStoreFactorySpi;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.Transaction;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureStore;
import org.geotools.feature.DefaultFeatureCollection;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.opengis.feature.simple.SimpleFeatureType;

@Description(OMSSHAPEFILEFEATUREWRITER_DESCRIPTION)
@Author(name = OMSSHAPEFILEFEATUREWRITER_AUTHORNAMES, contact = OMSSHAPEFILEFEATUREWRITER_AUTHORCONTACTS)
@Keywords(OMSSHAPEFILEFEATUREWRITER_KEYWORDS)
@Label(OMSSHAPEFILEFEATUREWRITER_LABEL)
@Name(OMSSHAPEFILEFEATUREWRITER_NAME)
@Status(OMSSHAPEFILEFEATUREWRITER_STATUS)
@License(OMSSHAPEFILEFEATUREWRITER_LICENSE)
@UI(OMSSHAPEFILEFEATUREWRITER_UI)
public class OmsShapefileFeatureWriter extends JGTModel {

    @Description(OMSSHAPEFILEFEATUREWRITER_geodata_DESCRIPTION)
    @In
    public SimpleFeatureCollection geodata = null;

    @Description(OMSSHAPEFILEFEATUREWRITER_file_DESCRIPTION)
    @UI(JGTConstants.FILEOUT_UI_HINT)
    @In
    public String file = null;

    @Description(OMSSHAPEFILEFEATUREWRITER_doIndex_DESCRIPTION)
    @In
    public boolean doIndex = true;

    @Description(OMSSHAPEFILEFEATUREWRITER_pType_DESCRIPTION)
    @In
    public SimpleFeatureType pType = null;

    private boolean hasWritten = false;

    @Execute
    public void writeFeatureCollection() throws IOException {
        if (!concatOr(!hasWritten, doReset)) {
            return;
        }

        pm.beginTask("Writing features to shapefile...", -1);

        if (!file.endsWith(".shp")) {
            file = file + ".shp";
        }
        if (geodata != null && geodata.size() != 0) {
            pType = geodata.getSchema();
        }
        File shapeFile = new File(file);
        FileDataStoreFactorySpi factory = FileDataStoreFinder.getDataStoreFactory("shp");
        Map map = Collections.singletonMap("url", shapeFile.toURI().toURL());
        DataStore newDataStore = factory.createNewDataStore(map);
        newDataStore.createSchema(pType);

        Transaction transaction = new DefaultTransaction("create");
        String typeName = newDataStore.getTypeNames()[0];
        SimpleFeatureStore featureStore = (SimpleFeatureStore) newDataStore.getFeatureSource(typeName);

        featureStore.setTransaction(transaction);
        try {
            if (geodata == null) {
                featureStore.addFeatures(new DefaultFeatureCollection());
            } else {
                featureStore.addFeatures(geodata);
            }
            transaction.commit();
        } catch (Exception problem) {
            transaction.rollback();
            throw new IOException(problem.getLocalizedMessage());
        } finally {
            transaction.close();
            pm.done();
        }

        hasWritten = true;
    }

    public static void writeShapefile( String path, SimpleFeatureCollection featureCollection, IJGTProgressMonitor pm )
            throws IOException {
        OmsShapefileFeatureWriter writer = new OmsShapefileFeatureWriter();
        if (pm != null) {
            writer.pm = pm;
        }
        writer.file = path;
        writer.geodata = featureCollection;
        writer.writeFeatureCollection();
    }

    public static void writeEmptyShapefile( String path, SimpleFeatureType schema, IJGTProgressMonitor pm ) throws IOException {
        OmsShapefileFeatureWriter writer = new OmsShapefileFeatureWriter();
        if (pm != null) {
            writer.pm = pm;
        }
        writer.file = path;
        writer.pType = schema;
        writer.writeFeatureCollection();
    }

}
