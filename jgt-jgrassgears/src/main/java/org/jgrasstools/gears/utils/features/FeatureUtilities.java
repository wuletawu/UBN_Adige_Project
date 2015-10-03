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
package org.jgrasstools.gears.utils.features;

import java.awt.geom.AffineTransform;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.text.MessageFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import javax.media.jai.JAI;
import javax.media.jai.ParameterBlockJAI;
import javax.media.jai.RenderedOp;
import javax.media.jai.iterator.RandomIter;

import org.geotools.coverage.grid.GridCoordinates2D;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.data.DataUtilities;
import org.geotools.data.DefaultTransaction;
import org.geotools.data.Transaction;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.simple.SimpleFeatureStore;
import org.geotools.feature.DefaultFeatureCollection;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.gce.grassraster.JGrassConstants;
import org.geotools.geometry.DirectPosition2D;
import org.geotools.geometry.Envelope2D;
import org.jaitools.media.jai.vectorize.VectorizeDescriptor;
import org.jgrasstools.gears.libs.exceptions.ModelsIOException;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.utils.RegionMap;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.gears.utils.geometry.GeometryType;
import org.jgrasstools.gears.utils.geometry.GeometryUtilities;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.feature.type.AttributeDescriptor;
import org.opengis.metadata.spatial.PixelOrientation;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.CoordinateList;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.LinearRing;
import com.vividsolutions.jts.geom.MultiLineString;
import com.vividsolutions.jts.geom.MultiPoint;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.geom.prep.PreparedGeometry;
import com.vividsolutions.jts.geom.prep.PreparedGeometryFactory;
import com.vividsolutions.jts.geom.util.AffineTransformation;
import com.vividsolutions.jts.index.strtree.STRtree;

public class FeatureUtilities {

    /**
     * Order the geometries of a list to be all directed in the same direction
     * 
     * @param geometryList the list of geometries to be ordered
     * @param thresHold a scalar value that defines the max distance between two points to be the
     *        same
     * @return a list of ordered coordinates
     */
    @SuppressWarnings("unchecked")
    public static CoordinateList orderLineGeometries( List<Geometry> geometryList, double thresHold ) {
        /*
         * first search the feature that is one of the two external points
         */
        Geometry firstFeature = null;
        boolean foundFirst = true;
        boolean foundSecond = true;
        for( Geometry feature : geometryList ) {
            foundFirst = true;
            foundSecond = true;

            Coordinate[] coords = feature.getCoordinates();

            Coordinate first = coords[0];
            Coordinate last = coords[coords.length - 1];

            for( Geometry compareFeature : geometryList ) {
                if (compareFeature.equals(feature))
                    continue;
                Coordinate[] compareCoords = compareFeature.getCoordinates();

                Coordinate comparefirst = compareCoords[0];
                Coordinate comparelast = compareCoords[compareCoords.length - 1];

                /*
                 * check if the next point is far away
                 */
                if (first.distance(comparefirst) < thresHold || first.distance(comparelast) < thresHold) {
                    foundFirst = false;
                }
                if (last.distance(comparefirst) < thresHold || last.distance(comparelast) < thresHold) {
                    foundSecond = false;
                }

            }
            if (foundFirst || foundSecond) {
                firstFeature = feature;
                break;
            }

        }
        if (firstFeature == null) {
            throw new RuntimeException();
        }
        CoordinateList coordinateList = new CoordinateList();
        Coordinate[] coords = firstFeature.getCoordinates();
        if (foundSecond) {
            for( int i = 0; i < coords.length; i++ ) {
                coordinateList.add(coords[coords.length - i - 1]);
            }
        } else {
            for( int i = 0; i < coords.length; i++ ) {
                coordinateList.add(coords[i]);
            }
        }

        // if (foundFirst) {
        // addCoordsInProperDirection(foundFirst, coordinateList, coords, true,
        // 0);
        // }else{
        // addCoordsInProperDirection(foundSecond, coordinateList, coords, true,
        // 0);
        // }

        geometryList.remove(firstFeature);

        Coordinate currentCoordinate = coordinateList.getCoordinate(coordinateList.size() - 1);
        while( geometryList.size() != 0 ) {

            for( int j = 0; j < geometryList.size(); j++ ) {
                System.out.println(j);
                Geometry compareGeom = geometryList.get(j);
                Coordinate[] compareCoords = compareGeom.getCoordinates();

                Coordinate comparefirst = compareCoords[0];
                Coordinate comparelast = compareCoords[compareCoords.length - 1];

                // System.out.println(j + " "
                // + currentCoordinate.distance(comparefirst) + " "
                // + currentCoordinate.distance(comparelast));

                /*
                 * check if the next point is far away
                 */
                if (currentCoordinate.distance(comparefirst) < thresHold) {
                    for( int i = 0; i < compareCoords.length; i++ ) {
                        coordinateList.add(compareCoords[i]);
                    }
                    currentCoordinate = new Coordinate(comparelast);
                    geometryList.remove(compareGeom);
                    break;
                } else if (currentCoordinate.distance(comparelast) < thresHold) {
                    for( int i = 0; i < compareCoords.length; i++ ) {
                        coordinateList.add(compareCoords[compareCoords.length - i - 1]);
                    }
                    currentCoordinate = new Coordinate(comparefirst);
                    geometryList.remove(compareGeom);
                    break;
                }

            }

        }

        return coordinateList;
    }

    /**
     * Create a featurecollection from a vector of features
     * 
     * @param features - the vectore of features
     * @return the created featurecollection
     */
    public static SimpleFeatureCollection createFeatureCollection( SimpleFeature... features ) {
        DefaultFeatureCollection fcollection = new DefaultFeatureCollection();

        for( SimpleFeature feature : features ) {
            fcollection.add(feature);
        }
        return fcollection;
    }

    /**
     * <p>
     * Convert a csv file to a FeatureCollection. 
     * <b>This for now supports only point geometries</b>.<br>
     * For different crs it also performs coor transformation.
     * </p>
     * <p>
     * <b>NOTE: this doesn't support date attributes</b>
     * </p>
     * 
     * @param csvFile the csv file.
     * @param crs the crs to use.
     * @param fieldsAndTypes the {@link Map} of filed names and {@link JGrassConstants#CSVTYPESARRAY types}.
     * @param pm progress monitor.
     * @param separatorthe separator to use, if null, comma is used.
     * @return the created {@link FeatureCollection}
     * @throws Exception
     */
    @SuppressWarnings("nls")
    public static SimpleFeatureCollection csvFileToFeatureCollection( File csvFile, CoordinateReferenceSystem crs,
            LinkedHashMap<String, Integer> fieldsAndTypesIndex, String separator, IJGTProgressMonitor pm ) throws Exception {
        GeometryFactory gf = new GeometryFactory();
        Map<String, Class< ? >> typesMap = JGrassConstants.CSVTYPESCLASSESMAP;
        String[] typesArray = JGrassConstants.CSVTYPESARRAY;

        if (separator == null) {
            separator = ",";
        }

        SimpleFeatureTypeBuilder b = new SimpleFeatureTypeBuilder();
        b.setName("csvimport");
        b.setCRS(crs);
        b.add("the_geom", Point.class);

        int xIndex = -1;
        int yIndex = -1;
        Set<String> fieldNames = fieldsAndTypesIndex.keySet();
        String[] fieldNamesArray = (String[]) fieldNames.toArray(new String[fieldNames.size()]);
        for( int i = 0; i < fieldNamesArray.length; i++ ) {
            String fieldName = fieldNamesArray[i];
            Integer typeIndex = fieldsAndTypesIndex.get(fieldName);

            if (typeIndex == 0) {
                xIndex = i;
            } else if (typeIndex == 1) {
                yIndex = i;
            } else {
                Class< ? > class1 = typesMap.get(typesArray[typeIndex]);
                b.add(fieldName, class1);
            }
        }
        SimpleFeatureType featureType = b.buildFeatureType();

        DefaultFeatureCollection newCollection = new DefaultFeatureCollection();
        Collection<Integer> orderedTypeIndexes = fieldsAndTypesIndex.values();
        Integer[] orderedTypeIndexesArray = (Integer[]) orderedTypeIndexes.toArray(new Integer[orderedTypeIndexes.size()]);

        BufferedReader bR = null;
        try {
            bR = new BufferedReader(new FileReader(csvFile));
            String line = null;
            int featureId = 0;
            pm.beginTask("Importing raw data", -1);
            while( (line = bR.readLine()) != null ) {
                pm.worked(1);
                if (line.startsWith("#")) {
                    continue;
                }

                SimpleFeatureBuilder builder = new SimpleFeatureBuilder(featureType);
                Object[] values = new Object[fieldNames.size() - 1];

                String[] lineSplit = line.split(separator);
                double x = Double.parseDouble(lineSplit[xIndex]);
                double y = Double.parseDouble(lineSplit[yIndex]);
                Point point = gf.createPoint(new Coordinate(x, y));
                values[0] = point;

                int objIndex = 1;
                for( int i = 0; i < lineSplit.length; i++ ) {
                    if (i == xIndex || i == yIndex) {
                        continue;
                    }

                    String value = lineSplit[i];
                    int typeIndex = orderedTypeIndexesArray[i];
                    String typeName = typesArray[typeIndex];
                    if (typeName.equals(typesArray[3])) {
                        values[objIndex] = value;
                    } else if (typeName.equals(typesArray[4])) {
                        values[objIndex] = new Double(value);
                    } else if (typeName.equals(typesArray[5])) {
                        values[objIndex] = new Integer(value);
                    } else {
                        throw new IllegalArgumentException("An undefined value type was found");
                    }
                    objIndex++;
                }
                builder.addAll(values);

                SimpleFeature feature = builder.buildFeature(featureType.getTypeName() + "." + featureId);
                featureId++;
                newCollection.add(feature);
            }
        } finally {
            if (bR != null)
                bR.close();
        }
        pm.done();

        return newCollection;
    }

    /**
     * <p>
     * The easy way to create a shapefile from attributes and geometries
     * </p>
     * <p>
     * <b>NOTE: this doesn't support date attributes</b>
     * </p>
     * 
     * @param shapeFilePath the shapefile name
     * @param crs the destination crs
     * @param fet the featurecollection
     * @throws IOException 
     */
    public static boolean collectionToShapeFile( String shapeFilePath, CoordinateReferenceSystem crs, SimpleFeatureCollection fet )
            throws IOException {

        // Create the file you want to write to
        File file = null;
        if (shapeFilePath.toLowerCase().endsWith(".shp")) { //$NON-NLS-1$
            file = new File(shapeFilePath);
        } else {
            file = new File(shapeFilePath + ".shp"); //$NON-NLS-1$
        }

        ShapefileDataStoreFactory factory = new ShapefileDataStoreFactory();

        Map<String, Serializable> create = new HashMap<String, Serializable>();
        create.put("url", file.toURI().toURL());
        ShapefileDataStore newDataStore = (ShapefileDataStore) factory.createNewDataStore(create);

        newDataStore.createSchema(fet.getSchema());
        if (crs != null)
            newDataStore.forceSchemaCRS(crs);
        Transaction transaction = new DefaultTransaction();
        SimpleFeatureStore featureStore = (SimpleFeatureStore) newDataStore.getFeatureSource();
        featureStore.setTransaction(transaction);
        try {
            featureStore.addFeatures(fet);
            transaction.commit();
        } catch (Exception problem) {
            problem.printStackTrace();
            transaction.rollback();
        } finally {
            transaction.close();
        }
        return true;
    }
    /**
     * @param name the shapefile name
     * @param fieldsSpec to create other fields you can use a string like : <br>
     *        "geom:MultiLineString,FieldName:java.lang.Integer" <br>
     *        field name can not be over 10 characters use a ',' between each field <br>
     *        field types can be : java.lang.Integer, java.lang.Long, // java.lang.Double,
     *        java.lang.String or java.util.Date
     * @return
     * @throws Exception 
     */
    @SuppressWarnings("nls")
    public static ShapefileDataStore createShapeFileDatastore( String name, String fieldsSpec, CoordinateReferenceSystem crs )
            throws Exception {
        // Create the file you want to write to
        File file = null;
        if (name.toLowerCase().endsWith(".shp")) {
            file = new File(name);
        } else {
            file = new File(name + ".shp");
        }

        ShapefileDataStoreFactory factory = new ShapefileDataStoreFactory();
        Map<String, Serializable> create = new HashMap<String, Serializable>();
        create.put("url", file.toURI().toURL());
        ShapefileDataStore myData = (ShapefileDataStore) factory.createNewDataStore(create);

        // Tell this shapefile what type of data it will store
        // Shapefile handle only : Point, MultiPoint, MultiLineString,
        // MultiPolygon
        SimpleFeatureType featureType = DataUtilities.createType(name, fieldsSpec);

        // Create the Shapefile (empty at this point)
        myData.createSchema(featureType);

        // Tell the DataStore what type of Coordinate Reference System (CRS)
        // to use
        myData.forceSchemaCRS(crs);

        return myData;

    }

    /**
     * Extracts features from a {@link FeatureCollection} into an {@link ArrayList}.
     * 
     * @param collection the feature collection.
     * @return the list with the features or an empty list if no features present.
     */
    public static List<SimpleFeature> featureCollectionToList( SimpleFeatureCollection collection ) {
        List<SimpleFeature> featuresList = new ArrayList<SimpleFeature>();
        if (collection == null) {
            return featuresList;
        }
        SimpleFeatureIterator featureIterator = collection.features();
        while( featureIterator.hasNext() ) {
            SimpleFeature feature = featureIterator.next();
            featuresList.add(feature);
        }
        featureIterator.close();
        return featuresList;
    }

    /**
     * Extracts features from a {@link FeatureCollection} into an {@link STRtree}.
     * 
     * @param collection the feature collection.
     * @return the tree containing the features.
     */
    public static STRtree featureCollectionToSTRtree( SimpleFeatureCollection collection ) {
        STRtree tree = new STRtree();
        SimpleFeatureIterator featureIterator = collection.features();
        while( featureIterator.hasNext() ) {
            SimpleFeature feature = featureIterator.next();
            Geometry geometry = (Geometry) feature.getDefaultGeometry();
            tree.insert(geometry.getEnvelopeInternal(), feature);
        }
        featureIterator.close();
        return tree;
    }

    /**
     * Extracts features from a {@link FeatureCollection} into an {@link ArrayList} of {@link FeatureMate}s.
     * 
     * @param collection the feature collection.
     * @return the list with the features or an empty list if no features present.
     */
    public static List<FeatureMate> featureCollectionToMatesList( SimpleFeatureCollection collection ) {
        List<FeatureMate> featuresList = new ArrayList<FeatureMate>();
        if (collection == null) {
            return featuresList;
        }
        SimpleFeatureIterator featureIterator = collection.features();
        while( featureIterator.hasNext() ) {
            SimpleFeature feature = featureIterator.next();
            featuresList.add(new FeatureMate(feature));
        }
        featureIterator.close();
        return featuresList;
    }

    /**
     * Extracts features from a {@link FeatureCollection} into an {@link ArrayList} of its geometries.
     * 
     * @param collection the feature collection.
     * @param doSubGeoms split the geometries in single geometries (ex. MultiLines in Lines).
     * @param userDataField if not <code>null</code>, the data in the field are put in the userData
     *                  field of the geometry.
     * @return the list with the geometries or an empty list if no features present.
     */
    public static List<Geometry> featureCollectionToGeometriesList( SimpleFeatureCollection collection, boolean doSubGeoms,
            String userDataField ) {
        List<Geometry> geometriesList = new ArrayList<Geometry>();
        if (collection == null) {
            return geometriesList;
        }
        SimpleFeatureIterator featureIterator = collection.features();
        while( featureIterator.hasNext() ) {
            SimpleFeature feature = featureIterator.next();
            Geometry geometry = (Geometry) feature.getDefaultGeometry();
            if (geometry != null)
                if (doSubGeoms) {
                    int numGeometries = geometry.getNumGeometries();
                    for( int i = 0; i < numGeometries; i++ ) {
                        Geometry geometryN = geometry.getGeometryN(i);
                        geometriesList.add(geometryN);
                        if (userDataField != null) {
                            Object attribute = feature.getAttribute(userDataField);
                            geometryN.setUserData(attribute);
                        }
                    }
                } else {
                    geometriesList.add(geometry);
                    if (userDataField != null) {
                        Object attribute = feature.getAttribute(userDataField);
                        geometry.setUserData(attribute);
                    }
                }
        }
        featureIterator.close();
        return geometriesList;
    }

    /**
     * Make a {@link SimpleFeatureCollection} from a set of {@link Geometry}.
     * 
     * <p>This is a fast utility and adds no attributes.</p>
     * 
     * @param crs The {@link CoordinateReferenceSystem}.
     * @param geometries the set of {@link Geometry} to add.
     * @return the features wrapping the geoms.
     */
    public static SimpleFeatureCollection featureCollectionFromGeometry( CoordinateReferenceSystem crs, Geometry... geometries ) {
        SimpleFeatureTypeBuilder b = new SimpleFeatureTypeBuilder();
        b.setName("simplegeom");
        b.setCRS(crs);

        DefaultFeatureCollection newCollection = new DefaultFeatureCollection();

        GeometryType geometryType = GeometryUtilities.getGeometryType(geometries[0]);
        switch( geometryType ) {
        case LINE:
            b.add("the_geom", LineString.class);
            break;
        case MULTILINE:
            b.add("the_geom", MultiLineString.class);
            break;
        case POINT:
            b.add("the_geom", Point.class);
            break;
        case MULTIPOINT:
            b.add("the_geom", MultiPoint.class);
            break;
        case POLYGON:
            b.add("the_geom", Polygon.class);
            break;
        case MULTIPOLYGON:
            b.add("the_geom", MultiPolygon.class);
            break;
        default:
            break;
        }

        Object userData = geometries[0].getUserData();
        int userDataSize = 1;
        if (userData != null) {
            if (userData instanceof String[]) {
                String[] string = (String[]) userData;
                userDataSize = string.length;
                for( int i = 0; i < userDataSize; i++ ) {
                    b.add("data" + i, String.class);
                }
            } else {
                b.add("userdata", userData.getClass());
            }
        }

        SimpleFeatureType type = b.buildFeatureType();
        SimpleFeatureBuilder builder = new SimpleFeatureBuilder(type);
        for( Geometry g : geometries ) {
            Object[] values;
            if (userData == null) {
                values = new Object[]{g};
            } else {
                Object tmpUserData = g.getUserData();
                if (tmpUserData instanceof String[]) {
                    String[] string = (String[]) tmpUserData;
                    values = new Object[userDataSize + 1];
                    values[0] = g;
                    for( int i = 0; i < string.length; i++ ) {
                        values[i + 1] = string[i];
                    }
                } else {
                    values = new Object[]{g, tmpUserData};
                }
            }
            builder.addAll(values);
            SimpleFeature feature = builder.buildFeature(null);
            newCollection.add(feature);
        }
        return newCollection;
    }

    /**
     * Getter for attributes of a feature.
     * 
     * <p>If the attribute is not found, checks are done in non
     * case sensitive mode.
     * 
     * @param feature the feature from which to get the attribute.
     * @param field the name of the field.
     * @return the attribute or null if none found.
     */
    public static Object getAttributeCaseChecked( SimpleFeature feature, String field ) {
        Object attribute = feature.getAttribute(field);
        if (attribute == null) {
            attribute = feature.getAttribute(field.toLowerCase());
            if (attribute != null)
                return attribute;
            attribute = feature.getAttribute(field.toUpperCase());
            if (attribute != null)
                return attribute;

            // alright, last try, search for it
            SimpleFeatureType featureType = feature.getFeatureType();
            field = findAttributeName(featureType, field);
            if (field != null) {
                return feature.getAttribute(field);
            }
        }
        return attribute;
    }

    /**
     * Find the name of an attribute, case insensitive.
     * 
     * @param featureType the feature type to check.
     * @param field the case insensitive field name.
     * @return the real name of the field, or <code>null</code>, if none found.
     */
    public static String findAttributeName( SimpleFeatureType featureType, String field ) {
        List<AttributeDescriptor> attributeDescriptors = featureType.getAttributeDescriptors();
        for( AttributeDescriptor attributeDescriptor : attributeDescriptors ) {
            String name = attributeDescriptor.getLocalName();
            if (name.toLowerCase().equals(field.toLowerCase())) {
                return name;
            }
        }
        return null;
    }

    /**
     * Create a {@link Polygon} from an {@link Envelope}.
     * 
     * @param envelope the envelope to convert.
     * @return the created polygon.
     */
    public static Polygon envelopeToPolygon( Envelope2D envelope ) {
        double w = envelope.getMinX();
        double e = envelope.getMaxX();
        double s = envelope.getMinY();
        double n = envelope.getMaxY();

        Coordinate[] coords = new Coordinate[5];
        coords[0] = new Coordinate(w, n);
        coords[1] = new Coordinate(e, n);
        coords[2] = new Coordinate(e, s);
        coords[3] = new Coordinate(w, s);
        coords[4] = new Coordinate(w, n);

        GeometryFactory gf = GeometryUtilities.gf();
        LinearRing linearRing = gf.createLinearRing(coords);
        Polygon polygon = gf.createPolygon(linearRing, null);
        return polygon;
    }

    /**
     * Create a {@link Polygon} from an {@link Envelope}.
     * 
     * @param envelope the envelope to convert.
     * @return the created polygon.
     */
    public static Polygon envelopeToPolygon( Envelope envelope ) {
        double w = envelope.getMinX();
        double e = envelope.getMaxX();
        double s = envelope.getMinY();
        double n = envelope.getMaxY();

        Coordinate[] coords = new Coordinate[5];
        coords[0] = new Coordinate(w, n);
        coords[1] = new Coordinate(e, n);
        coords[2] = new Coordinate(e, s);
        coords[3] = new Coordinate(w, s);
        coords[4] = new Coordinate(w, n);

        GeometryFactory gf = GeometryUtilities.gf();
        LinearRing linearRing = gf.createLinearRing(coords);
        Polygon polygon = gf.createPolygon(linearRing, null);
        return polygon;
    }

    /**
     * Extracts a list of polygons from the cell bounds of a given {@link GridCoverage2D coverage}.
     * 
     * <p><b>Note that the cells are added in a rows 
     * and cols order (for each row evaluate each column).</b></p> 
     * 
     * @param coverage the coverage to use.
     * @return the list of envelope geometries.
     */
    public static List<Polygon> gridcoverageToCellPolygons( GridCoverage2D coverage ) {
        RegionMap regionMap = CoverageUtilities.getRegionParamsFromGridCoverage(coverage);
        double west = regionMap.getWest();
        double north = regionMap.getNorth();
        double xres = regionMap.getXres();
        double yres = regionMap.getYres();
        int cols = regionMap.getCols();
        int rows = regionMap.getRows();

        List<Polygon> polygons = new ArrayList<Polygon>();
        for( int r = 0; r < rows; r++ ) {
            for( int c = 0; c < cols; c++ ) {
                double w = west + xres * c;
                double e = w + xres;
                double n = north - yres * r;
                double s = n - yres;

                Coordinate[] coords = new Coordinate[5];
                coords[0] = new Coordinate(w, n);
                coords[1] = new Coordinate(e, n);
                coords[2] = new Coordinate(e, s);
                coords[3] = new Coordinate(w, s);
                coords[4] = new Coordinate(w, n);

                GeometryFactory gf = GeometryUtilities.gf();
                LinearRing linearRing = gf.createLinearRing(coords);
                Polygon polygon = gf.createPolygon(linearRing, null);
                polygons.add(polygon);
            }
        }
        return polygons;
    }

    /**
     * Helper function to run the Vectorize operation with given parameters and
     * retrieve the vectors.
     * 
     * @param src the source {@link GridCoverage2D}.
     * @param args a {@code Map} of parameter names and values or <code>null</code>.
     * 
     * @return the generated vectors as JTS Polygons
     */
    @SuppressWarnings("unchecked")
    public static Collection<Polygon> doVectorize( GridCoverage2D src, Map<String, Object> args ) {
        if (args == null) {
            args = new HashMap<String, Object>();
        }

        ParameterBlockJAI pb = new ParameterBlockJAI("Vectorize");
        pb.setSource("source0", src.getRenderedImage());

        // Set any parameters that were passed in
        for( Entry<String, Object> e : args.entrySet() ) {
            pb.setParameter(e.getKey(), e.getValue());
        }

        // Get the desintation image: this is the unmodified source image data
        // plus a property for the generated vectors
        RenderedOp dest = JAI.create("Vectorize", pb);

        // Get the vectors
        Collection<Polygon> polygons = (Collection<Polygon>) dest.getProperty(VectorizeDescriptor.VECTOR_PROPERTY_NAME);

        RegionMap regionParams = CoverageUtilities.getRegionParamsFromGridCoverage(src);
        double xRes = regionParams.getXres();
        double yRes = regionParams.getYres();
        final AffineTransform mt2D = (AffineTransform) src.getGridGeometry().getGridToCRS2D(PixelOrientation.CENTER);
        final AffineTransformation jtsTransformation = new AffineTransformation(mt2D.getScaleX(), mt2D.getShearX(),
                mt2D.getTranslateX() - xRes / 2.0, mt2D.getShearY(), mt2D.getScaleY(), mt2D.getTranslateY() + yRes / 2.0);
        for( Polygon polygon : polygons ) {
            polygon.apply(jtsTransformation);
        }
        return polygons;
    }

    public static SimpleFeature toDummyFeature( Geometry geom, CoordinateReferenceSystem crs ) {
        SimpleFeatureTypeBuilder b = new SimpleFeatureTypeBuilder();
        b.setName("dummy");
        if (crs != null)
            b.setCRS(crs);
        b.add("the_geom", geom.getClass());
        SimpleFeatureType type = b.buildFeatureType();
        SimpleFeatureBuilder builder = new SimpleFeatureBuilder(type);
        Object[] values = new Object[]{geom};
        builder.addAll(values);
        SimpleFeature feature = builder.buildFeature(null);
        return feature;
    }

    /**
     * Extracts the positions and values of a polygon mapped on a raster(rasterization).
     * 
     * <p><b>Bound checks (polygon completely covering the raster) need to be
     * made before running this method.</b></p>
     * 
     * @param coverageIterator the raster from which to extract the values.
     * @param cols the cols of the raster.
     * @param rows the rows of the raster.
     * @param xRes the resolution of the raster.
     * @param gridGeometry the {@link GridGeometry2D} of the raster.
     * @param polygon the polygon to extract.
     * @param defaultValue a defaultvalue in case no raster is supplied.
     * @return the list of coordinates, with the z values set to the raster or default value.
     * @throws Exception
     */
    public static List<Coordinate> extractPolygonOnCoverage( RandomIter coverageIterator, int cols, int rows, double xRes,
            GridGeometry2D gridGeometry, Polygon polygon, double defaultValue ) throws Exception {
        GeometryFactory gf = GeometryUtilities.gf();
        PreparedGeometry preparedGeometry = PreparedGeometryFactory.prepare(polygon);
        double delta = xRes / 4.0;

        List<Coordinate> coordinatesList = new ArrayList<Coordinate>();

        for( int r = 0; r < rows; r++ ) {
            // do scan line to fill the polygon
            double[] westPos = gridGeometry.gridToWorld(new GridCoordinates2D(0, r)).getCoordinate();
            double[] eastPos = gridGeometry.gridToWorld(new GridCoordinates2D(cols - 1, r)).getCoordinate();
            Coordinate west = new Coordinate(westPos[0], westPos[1]);
            Coordinate east = new Coordinate(eastPos[0], eastPos[1]);
            LineString line = gf.createLineString(new Coordinate[]{west, east});
            if (preparedGeometry.intersects(line)) {
                Geometry internalLines = polygon.intersection(line);
                int lineNums = internalLines.getNumGeometries();
                for( int l = 0; l < lineNums; l++ ) {
                    Coordinate[] coords = internalLines.getGeometryN(l).getCoordinates();
                    if (coords.length == 2) {
                        for( int j = 0; j < coords.length; j = j + 2 ) {
                            Coordinate startC = new Coordinate(coords[j].x + delta, coords[j].y);
                            Coordinate endC = new Coordinate(coords[j + 1].x - delta, coords[j + 1].y);

                            DirectPosition2D startDP;
                            DirectPosition2D endDP;
                            if (startC.x < endC.x) {
                                startDP = new DirectPosition2D(startC.x, startC.x);
                                endDP = new DirectPosition2D(endC.x, endC.x);
                            } else {
                                startDP = new DirectPosition2D(endC.x, endC.x);
                                endDP = new DirectPosition2D(startC.x, startC.x);
                            }
                            GridCoordinates2D startGridCoord = gridGeometry.worldToGrid(startDP);
                            GridCoordinates2D endGridCoord = gridGeometry.worldToGrid(endDP);

                            /*
                             * the part in between has to be filled
                             */
                            for( int k = startGridCoord.x; k <= endGridCoord.x; k++ ) {
                                double v = defaultValue;
                                if (coverageIterator != null) {
                                    v = coverageIterator.getSampleDouble(k, r, 0);
                                }
                                double[] xy = gridGeometry.gridToWorld(new GridCoordinates2D(k, r)).getCoordinate();
                                Coordinate coordinate = new Coordinate(xy[0], xy[1], v);
                                coordinatesList.add(coordinate);
                            }
                        }
                    } else {
                        if (coords.length == 1) {
                            throw new ModelsIOException(
                                    MessageFormat.format("Found a cusp in: {0}/{1}", coords[0].x, coords[0].y),
                                    "FeatureUtilities");
                        } else {
                            throw new ModelsIOException(MessageFormat.format(
                                    "Found intersection with more than 2 points in: {0}/{1}", coords[0].x, coords[0].y),
                                    "FeatureUtilities");
                        }
                    }
                }

            }
        }

        return coordinatesList;
    }

    /**
     * Calculate the avg of a value in a list of {@link SimpleFeature}s.
     * 
     * <p>Empty records are ignored.
     * 
     * @param features the features.
     * @param field the field to consider. 
     * @return the avg.
     */
    public static double avg( List<SimpleFeature> features, String field ) {
        double sum = 0;
        int count = 0;
        for( SimpleFeature feature : features ) {
            Object attribute = feature.getAttribute(field);
            if (attribute instanceof Number) {
                sum = sum + ((Number) attribute).doubleValue();
                count++;
            }
        }
        double avg = sum / count;
        return avg;
    }

    /**
     * Calculate the histogram of a list of {@link SimpleFeature}s.
     * 
     * @param features the list of features.
     * @param field the field to consider. 
     * @param bins the number of bins.
     * @return the histogram as matrix of rows num like bins and 3 columns for [binCenter, count, cumulated-normalize-count].
     */
    public static double[][] histogram( List<SimpleFeature> features, String field, int bins ) {
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;

        for( SimpleFeature feature : features ) {
            Object attribute = feature.getAttribute(field);
            if (attribute instanceof Number) {
                double value = ((Number) attribute).doubleValue();
                min = Math.min(min, value);
                max = Math.max(max, value);
            }
        }

        double range = max - min;
        double step = range / bins;
        double[][] histogram = new double[bins][3];
        for( int i = 0; i < histogram.length; i++ ) {
            histogram[i][0] = min + step * (i + 1);
        }

        for( SimpleFeature feature : features ) {
            Object attribute = feature.getAttribute(field);
            if (attribute instanceof Number) {
                double value = ((Number) attribute).doubleValue();
                for( int j = 0; j < histogram.length; j++ ) {
                    if (value <= histogram[j][0]) {
                        histogram[j][1] = histogram[j][1] + 1;
                        break;
                    }
                }
            }
        }

        double cumulatedMax = 0;
        for( int i = 0; i < histogram.length; i++ ) {
            if (i == 0) {
                histogram[i][2] = histogram[i][1];
            } else {
                histogram[i][2] = (histogram[i - 1][2] + histogram[i][1]);
            }
            cumulatedMax = histogram[i][2];
        }

        for( int i = 0; i < histogram.length; i++ ) {
            histogram[i][2] = histogram[i][2] / cumulatedMax;
            // and move the bin markers to their centers
            histogram[i][0] = histogram[i][0] - step / 2.0;
        }

        return histogram;
    }

}
