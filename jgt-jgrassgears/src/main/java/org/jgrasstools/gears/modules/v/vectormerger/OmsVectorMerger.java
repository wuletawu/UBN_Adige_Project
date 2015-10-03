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
package org.jgrasstools.gears.modules.v.vectormerger;

import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORMERGER_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORMERGER_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORMERGER_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORMERGER_DOCUMENTATION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORMERGER_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORMERGER_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORMERGER_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORMERGER_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORMERGER_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORMERGER_inVectors_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORMERGER_outVector_DESCRIPTION;

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

import org.geotools.data.DataUtilities;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.feature.DefaultFeatureCollection;
import org.jgrasstools.gears.libs.exceptions.ModelsIllegalargumentException;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.utils.features.FeatureExtender;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

@Description(OMSVECTORMERGER_DESCRIPTION)
@Documentation(OMSVECTORMERGER_DOCUMENTATION)
@Author(name = OMSVECTORMERGER_AUTHORNAMES, contact = OMSVECTORMERGER_AUTHORCONTACTS)
@Keywords(OMSVECTORMERGER_KEYWORDS)
@Label(OMSVECTORMERGER_LABEL)
@Name(OMSVECTORMERGER_NAME)
@Status(OMSVECTORMERGER_STATUS)
@License(OMSVECTORMERGER_LICENSE)
public class OmsVectorMerger extends JGTModel {

    @Description(OMSVECTORMERGER_inVectors_DESCRIPTION)
    @In
    public List<SimpleFeatureCollection> inVectors;

    @Description(OMSVECTORMERGER_outVector_DESCRIPTION)
    @Out
    public SimpleFeatureCollection outVector;

    @Execute
    public void process() throws Exception {
        checkNull(inVectors);

        SimpleFeatureType firstType = null;

        FeatureExtender fEx = null;

        pm.beginTask("Merging features...", inVectors.size());
        try {
            outVector = new DefaultFeatureCollection();
            for( SimpleFeatureCollection featureCollection : inVectors ) {
                if (firstType == null) {
                    firstType = featureCollection.getSchema();
                    fEx = new FeatureExtender(firstType, new String[0], new Class< ? >[0]);
                } else {
                    SimpleFeatureType schema = featureCollection.getSchema();
                    int compare = DataUtilities.compare(firstType, schema);
                    if (compare != 0) {
                        throw new ModelsIllegalargumentException("Merging is done only on same feature types.", this, pm);
                    }
                }
                SimpleFeatureIterator featureIterator = featureCollection.features();
                while( featureIterator.hasNext() ) {
                    SimpleFeature f = featureIterator.next();

                    SimpleFeature extendFeature = fEx.extendFeature(f, new Object[0]);

                    ((DefaultFeatureCollection) outVector).add(extendFeature);
                }
                pm.worked(1);
            }
        } finally {
            pm.done();
        }
    }
}
