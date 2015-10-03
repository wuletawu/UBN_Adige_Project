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

import static org.jgrasstools.gears.i18n.GearsMessages.OMSCOVERAGEVIEWER_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCOVERAGEVIEWER_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCOVERAGEVIEWER_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCOVERAGEVIEWER_DOCUMENTATION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCOVERAGEVIEWER_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCOVERAGEVIEWER_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCOVERAGEVIEWER_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCOVERAGEVIEWER_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCOVERAGEVIEWER_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCOVERAGEVIEWER_UI;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSCOVERAGEVIEWER_raster_DESCRIPTION;

import java.awt.Color;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.RenderedImage;

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
import org.geotools.factory.CommonFactoryFinder;
import org.geotools.map.DefaultMapContext;
import org.geotools.map.MapContext;
import org.geotools.styling.ColorMap;
import org.geotools.styling.ColorMapEntry;
import org.geotools.styling.RasterSymbolizer;
import org.geotools.styling.SLD;
import org.geotools.styling.Style;
import org.geotools.styling.StyleBuilder;
import org.geotools.styling.StyleFactory;
import org.geotools.swing.JMapFrame;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.opengis.filter.expression.Expression;

@Description(OMSCOVERAGEVIEWER_DESCRIPTION)
@Documentation(OMSCOVERAGEVIEWER_DOCUMENTATION)
@Author(name = OMSCOVERAGEVIEWER_AUTHORNAMES, contact = OMSCOVERAGEVIEWER_AUTHORCONTACTS)
@Keywords(OMSCOVERAGEVIEWER_KEYWORDS)
@Label(OMSCOVERAGEVIEWER_LABEL)
@Name(OMSCOVERAGEVIEWER_NAME)
@Status(OMSCOVERAGEVIEWER_STATUS)
@License(OMSCOVERAGEVIEWER_LICENSE)
@UI(OMSCOVERAGEVIEWER_UI)
public class OmsCoverageViewer extends JGTModel {

    @Description(OMSCOVERAGEVIEWER_raster_DESCRIPTION)
    @In
    public GridCoverage2D raster = null;

    private StyleFactory sf = CommonFactoryFinder.getStyleFactory(null);

    @Execute
    public void viewCoverage() throws Exception {
        // RasterSymbolizer sym = sf.getDefaultRasterSymbolizer();
        // Style rasterStyle = SLD.wrapSymbolizers(sym);

        StyleBuilder sB = new StyleBuilder(sf);
        RasterSymbolizer rasterSym = sf.createRasterSymbolizer();

        ColorMap colorMap = sf.createColorMap();

        RenderedImage renderedImage = raster.getRenderedImage();
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
        Color toColor = Color.red;
        Expression fromColorExpr = sB.colorExpression(new java.awt.Color(fromColor.getRed(), fromColor.getGreen(), fromColor
                .getBlue(), 255));
        Expression toColorExpr = sB.colorExpression(new java.awt.Color(toColor.getRed(), toColor.getGreen(), toColor.getBlue(),
                255));
        Expression fromExpr = sB.literalExpression(min);
        Expression toExpr = sB.literalExpression(max);

        ColorMapEntry entry = sf.createColorMapEntry();
        entry.setQuantity(fromExpr);
        entry.setColor(fromColorExpr);
        colorMap.addColorMapEntry(entry);

        entry = sf.createColorMapEntry();
        entry.setQuantity(toExpr);
        entry.setColor(toColorExpr);
        colorMap.addColorMapEntry(entry);

        rasterSym.setColorMap(colorMap);

        Style rasterStyle = SLD.wrapSymbolizers(rasterSym);

        // Set up a MapContext with the two layers
        final MapContext map = new DefaultMapContext();
        map.setTitle("Coverage Viewer");
        map.addLayer(raster, rasterStyle);

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

}
