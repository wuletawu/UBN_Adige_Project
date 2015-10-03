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
package org.jgrasstools.gears.modules.v.vectorsimplifier;

import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORSIMPLIFIER_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORSIMPLIFIER_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORSIMPLIFIER_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORSIMPLIFIER_DOCUMENTATION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORSIMPLIFIER_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORSIMPLIFIER_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORSIMPLIFIER_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORSIMPLIFIER_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORSIMPLIFIER_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORSIMPLIFIER_inVector_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORSIMPLIFIER_outVector_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORSIMPLIFIER_pTolerance_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSVECTORSIMPLIFIER_pType_DESCRIPTION;

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

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.DefaultFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.utils.features.FeatureGeometrySubstitutor;
import org.jgrasstools.gears.utils.geometry.GeometryUtilities;
import org.opengis.feature.simple.SimpleFeature;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryCollection;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.simplify.DouglasPeuckerSimplifier;
import com.vividsolutions.jts.simplify.TopologyPreservingSimplifier;

@Description(OMSVECTORSIMPLIFIER_DESCRIPTION)
@Documentation(OMSVECTORSIMPLIFIER_DOCUMENTATION)
@Author(name = OMSVECTORSIMPLIFIER_AUTHORNAMES, contact = OMSVECTORSIMPLIFIER_AUTHORCONTACTS)
@Keywords(OMSVECTORSIMPLIFIER_KEYWORDS)
@Label(OMSVECTORSIMPLIFIER_LABEL)
@Name(OMSVECTORSIMPLIFIER_NAME)
@Status(OMSVECTORSIMPLIFIER_STATUS)
@License(OMSVECTORSIMPLIFIER_LICENSE)
public class OmsVectorSimplifier extends JGTModel {

    @Description(OMSVECTORSIMPLIFIER_inVector_DESCRIPTION)
    @In
    public SimpleFeatureCollection inVector;

    @Description(OMSVECTORSIMPLIFIER_pType_DESCRIPTION)
    @In
    public int pType = 0;

    @Description(OMSVECTORSIMPLIFIER_pTolerance_DESCRIPTION)
    @In
    public double pTolerance = 0.2;

    @Description(OMSVECTORSIMPLIFIER_outVector_DESCRIPTION)
    @Out
    public SimpleFeatureCollection outVector;

    private GeometryFactory gF = GeometryUtilities.gf();

    @Execute
    public void process() throws Exception {
        if (!concatOr(outVector == null, doReset)) {
            return;
        }
        FeatureIterator<SimpleFeature> inFeatureIterator = inVector.features();

        outVector = new DefaultFeatureCollection();

        FeatureGeometrySubstitutor fGS = new FeatureGeometrySubstitutor(inVector.getSchema());

        int size = inVector.size();
        pm.beginTask("Simplifing features...", size);
        while( inFeatureIterator.hasNext() ) {
            SimpleFeature feature = inFeatureIterator.next();

            Geometry geometry = (Geometry) feature.getDefaultGeometry();

            List<Geometry> geomList = new ArrayList<Geometry>();

            int numGeometries = geometry.getNumGeometries();
            for( int i = 0; i < numGeometries; i++ ) {
                Geometry geometryN = geometry.getGeometryN(i);
                switch( pType ) {
                case 0:
                    TopologyPreservingSimplifier tpSimplifier = new TopologyPreservingSimplifier(geometryN);
                    tpSimplifier.setDistanceTolerance(pTolerance);
                    Geometry tpsGeometry = tpSimplifier.getResultGeometry();
                    geomList.add(tpsGeometry);
                    break;
                case 1:
                    DouglasPeuckerSimplifier dpSimplifier = new DouglasPeuckerSimplifier(geometryN);
                    dpSimplifier.setDistanceTolerance(pTolerance);
                    Geometry dpsGeometry = dpSimplifier.getResultGeometry();
                    geomList.add(dpsGeometry);
                    break;
                default:
                }

            }

            Geometry newGeometry = null;
            if (geomList.size() == 1) {
                newGeometry = geomList.get(0);
            } else {
                Geometry[] geomArray = (Geometry[]) geomList.toArray(new Geometry[geomList.size()]);
                newGeometry = new GeometryCollection(geomArray, gF);
            }

            SimpleFeature newFeature = fGS.substituteGeometry(feature, newGeometry);

            ((DefaultFeatureCollection) outVector).add(newFeature);
            pm.worked(1);
        }
        pm.done();

        inFeatureIterator.close();

    }

}
