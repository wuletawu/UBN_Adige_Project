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
package org.jgrasstools.gears.ui;

import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_DOCUMENTATION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_UI;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_inRaster_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_inRasters_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_inSld_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_inVector_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSMAPSVIEWER_inVectors_DESCRIPTION;

import java.awt.Color;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.RenderedImage;
import java.io.File;

import javax.media.jai.iterator.RectIter;
import javax.media.jai.iterator.RectIterFactory;

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

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.factory.CommonFactoryFinder;
import org.geotools.factory.GeoTools;
import org.geotools.map.DefaultMapContext;
import org.geotools.map.MapContext;
import org.geotools.styling.ColorMap;
import org.geotools.styling.ColorMapEntry;
import org.geotools.styling.FeatureTypeStyle;
import org.geotools.styling.Fill;
import org.geotools.styling.Graphic;
import org.geotools.styling.Mark;
import org.geotools.styling.PointSymbolizer;
import org.geotools.styling.PolygonSymbolizer;
import org.geotools.styling.RasterSymbolizer;
import org.geotools.styling.Rule;
import org.geotools.styling.SLD;
import org.geotools.styling.SLDParser;
import org.geotools.styling.SLDTransformer;
import org.geotools.styling.Stroke;
import org.geotools.styling.Style;
import org.geotools.styling.StyleBuilder;
import org.geotools.styling.StyleFactory;
import org.geotools.styling.StyledLayerDescriptor;
import org.geotools.swing.JMapFrame;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.utils.geometry.GeometryType;
import org.jgrasstools.gears.utils.geometry.GeometryUtilities;
import org.opengis.feature.type.GeometryDescriptor;
import org.opengis.filter.FilterFactory;
import org.opengis.filter.expression.Expression;

@Description(OMSMAPSVIEWER_DESCRIPTION)
@Documentation(OMSMAPSVIEWER_DOCUMENTATION)
@Author(name = OMSMAPSVIEWER_AUTHORNAMES, contact = OMSMAPSVIEWER_AUTHORCONTACTS)
@Keywords(OMSMAPSVIEWER_KEYWORDS)
@Label(OMSMAPSVIEWER_LABEL)
@Name(OMSMAPSVIEWER_NAME)
@Status(OMSMAPSVIEWER_STATUS)
@License(OMSMAPSVIEWER_LICENSE)
@UI(OMSMAPSVIEWER_UI)
public class OmsMapsViewer extends JGTModel {

    @Description(OMSMAPSVIEWER_inRasters_DESCRIPTION)
    @In
    public GridCoverage2D[] inRasters = new GridCoverage2D[0];

    @Description(OMSMAPSVIEWER_inRaster_DESCRIPTION)
    @In
    public GridCoverage2D inRaster = null;

    @Description(OMSMAPSVIEWER_inVectors_DESCRIPTION)
    @In
    public SimpleFeatureCollection[] inVectors = new SimpleFeatureCollection[0];

    @Description(OMSMAPSVIEWER_inVector_DESCRIPTION)
    @In
    public SimpleFeatureCollection inVector = null;

    @Description(OMSMAPSVIEWER_inSld_DESCRIPTION)
    @In
    public String inSld = null;

    private StyleFactory sf = CommonFactoryFinder.getStyleFactory(GeoTools.getDefaultHints());
    private FilterFactory ff = CommonFactoryFinder.getFilterFactory(GeoTools.getDefaultHints());
    private StyleBuilder sb = new StyleBuilder(sf, ff);

    private Style namedStyle;

    @Execute
    public void displayMaps() throws Exception {
        final MapContext map = new DefaultMapContext();
        map.setTitle("Maps Viewer");

        RasterSymbolizer rasterSym = sf.createRasterSymbolizer();

        if (inRaster != null) {
            inRasters = new GridCoverage2D[]{inRaster};
        }
        addCoverages(map, sb, rasterSym);

        if (inVector != null) {
            inVectors = new SimpleFeatureCollection[]{inVector};
            // does it have style
            if (inSld != null) {
                File sldFile = new File(inSld);
                if (sldFile.exists()) {
                    SLDParser stylereader = new SLDParser(sf, sldFile);
                    StyledLayerDescriptor sld = stylereader.parseSLD();

                    namedStyle = SLD.defaultStyle(sld);
                    SLDTransformer aTransformer = new SLDTransformer();
                    aTransformer.setIndentation(4);
                    String xml = ""; //$NON-NLS-1$
                    xml = aTransformer.transform(sld);
                    System.out.println(xml);
                }
            }

        }
        addFeatureCollections(map);

        // Create a JMapFrame with a menu to choose the display style for the
        final JMapFrame frame = new JMapFrame(map);
        frame.setSize(800, 600);
        frame.enableStatusBar(true);
        frame.enableTool(JMapFrame.Tool.ZOOM, JMapFrame.Tool.PAN, JMapFrame.Tool.RESET);
        frame.enableToolBar(true);
        frame.setVisible(true);
        frame.addWindowListener(new WindowAdapter(){
            public void windowClosing( WindowEvent e ) {
                frame.setVisible(false);
            }
        });

        while( frame.isVisible() ) {
            Thread.sleep(300);
        }
    }

    private void addFeatureCollections( MapContext map ) {
        for( SimpleFeatureCollection fc : inVectors ) {
            GeometryDescriptor geometryDescriptor = fc.getSchema().getGeometryDescriptor();
            GeometryType type = GeometryUtilities.getGeometryType(geometryDescriptor.getType());

            switch( type ) {
            case MULTIPOLYGON:
            case POLYGON:

                Stroke polygonStroke = sf.createStroke(ff.literal(Color.BLACK), ff.literal(1));
                Fill polygonFill = sf.createFill(ff.literal(Color.RED), ff.literal(0.5));

                Rule polygonRule = sf.createRule();
                PolygonSymbolizer polygonSymbolizer = sf.createPolygonSymbolizer(polygonStroke, polygonFill, null);
                polygonRule.symbolizers().add(polygonSymbolizer);

                FeatureTypeStyle polygonFeatureTypeStyle = sf.createFeatureTypeStyle();
                polygonFeatureTypeStyle.rules().add(polygonRule);

                namedStyle = sf.createStyle();
                namedStyle.featureTypeStyles().add(polygonFeatureTypeStyle);
                namedStyle.setName("polygons");

                break;
            case MULTIPOINT:
            case POINT:
                if (namedStyle == null) {
                    Mark circleMark = sf.getCircleMark();
                    Fill fill = sf.createFill(ff.literal(Color.RED));
                    circleMark.setFill(fill);
                    // circleMark.setStroke(null);

                    Graphic gr = sf.createDefaultGraphic();
                    gr.graphicalSymbols().clear();
                    gr.graphicalSymbols().add(circleMark);
                    Expression size = ff.literal(6);
                    gr.setSize(size);

                    Rule pointRule = sf.createRule();
                    PointSymbolizer pointSymbolizer = sf.createPointSymbolizer(gr, null);
                    pointRule.symbolizers().add(pointSymbolizer);

                    FeatureTypeStyle pointsFeatureTypeStyle = sf.createFeatureTypeStyle();
                    pointsFeatureTypeStyle.rules().add(pointRule);

                    namedStyle = sf.createStyle();
                    namedStyle.featureTypeStyles().add(pointsFeatureTypeStyle);
                    namedStyle.setName("points");
                }
                break;
            case MULTILINE:
            case LINE:

                break;

            default:
                break;
            }

            map.addLayer(fc, namedStyle);

        }

    }

    private void addCoverages( final MapContext map, StyleBuilder sB, RasterSymbolizer rasterSym ) {
        ColorMap colorMap = sf.createColorMap();

        for( GridCoverage2D coverage : inRasters ) {
            RenderedImage renderedImage = coverage.getRenderedImage();
            double max = Double.NEGATIVE_INFINITY;
            double min = Double.POSITIVE_INFINITY;
            RectIter iter = RectIterFactory.create(renderedImage, null);
            do {
                do {
                    double value = iter.getSampleDouble();
                    if (value > max) {
                        max = value;
                    }
                    if (value < min) {
                        min = value;
                    }
                } while( !iter.nextPixelDone() );
                iter.startPixels();
            } while( !iter.nextLineDone() );

            // red to blue
            Color fromColor = Color.blue;
            Color midColor = Color.green;
            Color toColor = Color.red;
            Expression fromColorExpr = sB.colorExpression(new java.awt.Color(fromColor.getRed(), fromColor.getGreen(), fromColor
                    .getBlue(), 255));
            Expression midColorExpr = sB.colorExpression(new java.awt.Color(midColor.getRed(), midColor.getGreen(), midColor
                    .getBlue(), 255));
            Expression toColorExpr = sB.colorExpression(new java.awt.Color(toColor.getRed(), toColor.getGreen(), toColor
                    .getBlue(), 255));
            Expression fromExpr = sB.literalExpression(min);
            Expression midExpr = sB.literalExpression(min + (max - min) / 2);
            Expression toExpr = sB.literalExpression(max);

            ColorMapEntry entry = sf.createColorMapEntry();
            entry.setQuantity(fromExpr);
            entry.setColor(fromColorExpr);
            colorMap.addColorMapEntry(entry);

            entry = sf.createColorMapEntry();
            entry.setQuantity(midExpr);
            entry.setColor(midColorExpr);
            colorMap.addColorMapEntry(entry);

            entry = sf.createColorMapEntry();
            entry.setQuantity(toExpr);
            entry.setColor(toColorExpr);
            colorMap.addColorMapEntry(entry);

            rasterSym.setColorMap(colorMap);

            Style rasterStyle = SLD.wrapSymbolizers(rasterSym);

            map.addLayer(coverage, rasterStyle);
        }
    }

    public static synchronized void displayRasterAndFeatures( final GridCoverage2D raster,
            final SimpleFeatureCollection... vectors ) throws Exception {

        new Thread(){
            @Override
            public void run() {
                OmsMapsViewer viewer = new OmsMapsViewer();
                if (raster != null) {
                    viewer.inRasters = new GridCoverage2D[]{raster};
                }
                if (vectors != null) {
                    viewer.inVectors = vectors;
                }
                try {
                    viewer.displayMaps();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }.start();

    }

    // @SuppressWarnings("nls")
    // public static void main( String[] args ) throws Exception {
    // GridCoverage2D coverage = OmsRasterReader.readRaster("/home/moovida/TMP/byumba_basins.asc");
    // SimpleFeatureCollection shapefile = OmsShapefileFeatureReader
    // .readShapefile("/home/moovida/TMP/byumba_extrbasins.shp");
    // displayRasterAndFeatures(coverage, shapefile);
    // }

}
