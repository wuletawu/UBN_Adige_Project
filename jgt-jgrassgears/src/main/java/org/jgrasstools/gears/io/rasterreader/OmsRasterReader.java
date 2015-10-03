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
package org.jgrasstools.gears.io.rasterreader;

import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_doLegacyGrass_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_fileNovalue_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_file_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_geodataNovalue_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_outRaster_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_pCols_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_pEast_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_pNorth_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_pRows_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_pSouth_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_pWest_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_pXres_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSRASTERREADER_pYres_DESCRIPTION;
import static org.jgrasstools.gears.libs.modules.JGTConstants.AIG;
import static org.jgrasstools.gears.libs.modules.JGTConstants.ESRIGRID;
import static org.jgrasstools.gears.libs.modules.JGTConstants.GEOTIF;
import static org.jgrasstools.gears.libs.modules.JGTConstants.GEOTIFF;
import static org.jgrasstools.gears.libs.modules.JGTConstants.GRASS;
import static org.jgrasstools.gears.libs.modules.JGTConstants.JPEG;
import static org.jgrasstools.gears.libs.modules.JGTConstants.JPG;
import static org.jgrasstools.gears.libs.modules.JGTConstants.PNG;
import static org.jgrasstools.gears.libs.modules.JGTConstants.doubleNovalue;
import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.COLS;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.EAST;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.NORTH;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.ROWS;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.SOUTH;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.WEST;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.XRES;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.YRES;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.buildCoverage;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.createDoubleWritableRaster;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.createGridGeometryGeneralParameter;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.getRegionParamsFromGridCoverage;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.gridGeometryFromRegionParams;
import static org.jgrasstools.gears.utils.coverage.CoverageUtilities.makeRegionParamsMap;

import java.awt.RenderingHints;
import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.HashMap;

import javax.media.jai.ImageLayout;
import javax.media.jai.JAI;
import javax.media.jai.iterator.RandomIter;
import javax.media.jai.iterator.RandomIterFactory;
import javax.media.jai.iterator.WritableRandomIter;

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

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.coverage.grid.ViewType;
import org.geotools.coverage.grid.io.AbstractGridCoverage2DReader;
import org.geotools.coverage.processing.Operations;
import org.geotools.coverageio.gdal.BaseGDALGridCoverage2DReader;
import org.geotools.coverageio.gdal.aig.AIGReader;
import org.geotools.factory.Hints;
import org.geotools.gce.arcgrid.ArcGridReader;
import org.geotools.gce.geotiff.GeoTiffReader;
import org.geotools.gce.grassraster.GrassCoverageReader;
import org.geotools.gce.grassraster.JGrassMapEnvironment;
import org.geotools.gce.grassraster.JGrassRegion;
import org.geotools.gce.grassraster.format.GrassCoverageFormat;
import org.geotools.gce.grassraster.format.GrassCoverageFormatFactory;
import org.geotools.gce.image.WorldImageReader;
import org.geotools.geometry.GeneralEnvelope;
import org.geotools.geometry.jts.ReferencedEnvelope;
import org.jgrasstools.gears.io.grasslegacy.OmsGrassLegacyReader;
import org.jgrasstools.gears.io.grasslegacy.utils.GrassLegacyUtilities;
import org.jgrasstools.gears.libs.exceptions.ModelsIllegalargumentException;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.gears.utils.files.FileUtilities;
import org.jgrasstools.gears.utils.math.NumericsUtilities;
import org.opengis.parameter.GeneralParameterValue;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import com.vividsolutions.jts.geom.Envelope;

@Description(OMSRASTERREADER_DESCRIPTION)
@Author(name = OMSRASTERREADER_AUTHORNAMES, contact = OMSRASTERREADER_AUTHORCONTACTS)
@Keywords(OMSRASTERREADER_KEYWORDS)
@Label(OMSRASTERREADER_LABEL)
@Name(OMSRASTERREADER_NAME)
@Status(OMSRASTERREADER_STATUS)
@License(OMSRASTERREADER_LICENSE)
public class OmsRasterReader extends JGTModel {

    @Description(OMSRASTERREADER_file_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String file = null;

    @Description(OMSRASTERREADER_fileNovalue_DESCRIPTION)
    @In
    public Double fileNovalue = -9999.0;

    @Description(OMSRASTERREADER_geodataNovalue_DESCRIPTION)
    @In
    public Double geodataNovalue = doubleNovalue;

    @Description(OMSRASTERREADER_pNorth_DESCRIPTION)
    @UI(JGTConstants.PROCESS_NORTH_UI_HINT)
    @In
    public Double pNorth = null;

    @Description(OMSRASTERREADER_pSouth_DESCRIPTION)
    @UI(JGTConstants.PROCESS_SOUTH_UI_HINT)
    @In
    public Double pSouth = null;

    @Description(OMSRASTERREADER_pWest_DESCRIPTION)
    @UI(JGTConstants.PROCESS_WEST_UI_HINT)
    @In
    public Double pWest = null;

    @Description(OMSRASTERREADER_pEast_DESCRIPTION)
    @UI(JGTConstants.PROCESS_EAST_UI_HINT)
    @In
    public Double pEast = null;

    @Description(OMSRASTERREADER_pXres_DESCRIPTION)
    @UI(JGTConstants.PROCESS_XRES_UI_HINT)
    @In
    public Double pXres = null;

    @Description(OMSRASTERREADER_pYres_DESCRIPTION)
    @UI(JGTConstants.PROCESS_YRES_UI_HINT)
    @In
    public Double pYres = null;

    @Description(OMSRASTERREADER_pRows_DESCRIPTION)
    @UI(JGTConstants.PROCESS_ROWS_UI_HINT)
    @In
    public Integer pRows = null;

    @Description(OMSRASTERREADER_pCols_DESCRIPTION)
    @UI(JGTConstants.PROCESS_COLS_UI_HINT)
    @In
    public Integer pCols = null;

    @Description(OMSRASTERREADER_doLegacyGrass_DESCRIPTION)
    @In
    public Boolean doLegacyGrass = false;

    @Description(OMSRASTERREADER_outRaster_DESCRIPTION)
    @Out
    public GridCoverage2D outRaster = null;

    /**
     * Flag to read only envelope (if true, the output geodata is null).
     */
    public boolean doEnvelope = false;

    /**
     * The original envelope of the coverage.
     */
    public GeneralEnvelope originalEnvelope;

    private GeneralParameterValue[] generalParameter = null;

    private double[] pBounds;

    private double[] pRes;
    private int[] pRowcol;

    private double internalFileNovalue = -9999.0;
    private double internalGeodataNovalue = doubleNovalue;

    private boolean doLegacyGrassAutoBoxed = false;

    @Execute
    public void process() throws Exception {
        if (!concatOr(outRaster == null, doReset)) {
            return;
        }

        if (doLegacyGrass != null && doLegacyGrass) {
            doLegacyGrassAutoBoxed = true;
        }

        if (fileNovalue != null) {
            internalFileNovalue = fileNovalue;
        }
        if (geodataNovalue != null) {
            internalGeodataNovalue = geodataNovalue;
        }

        if (hasBoundsRequest() && (!hasResolutionRequest() && !hasRowColsRequest())) {
            throw new RuntimeException("If bounds are requested, also a resolution or number of rows/cols has to be supplied.");
        }
        if (hasBoundsRequest()) {
            pBounds = new double[]{pNorth, pSouth, pWest, pEast};
        }
        if (hasResolutionRequest()) {
            pRes = new double[]{pXres, pYres};
        }
        if (hasRowColsRequest()) {
            pRowcol = new int[]{pRows, pCols};
        }

        String pType = null;
        // try to guess from the extension
        if (file.toLowerCase().endsWith(ESRIGRID)) {
            pType = ESRIGRID;
        } else if (file.toLowerCase().endsWith(AIG)) {
            pType = AIG;
        } else if (file.toLowerCase().endsWith(GEOTIFF) || file.toLowerCase().endsWith(GEOTIF)) {
            pType = GEOTIFF;
        } else if (file.toLowerCase().endsWith(JPEG) || file.toLowerCase().endsWith(JPG)) {
            pType = JPG;
        } else if (file.toLowerCase().endsWith(PNG)) {
            pType = PNG;
        } else if (CoverageUtilities.isGrass(file)) {
            pType = GRASS;
        } else
            throw new ModelsIllegalargumentException(
                    "Can't recognize the data format. Supported are: asc, tiff, jpg, png, grass.", this.getClass()
                            .getSimpleName(), pm);

        File mapFile = new File(file);
        try {
            pm.beginTask("Reading coverage: " + mapFile.getName(), IJGTProgressMonitor.UNKNOWN);

            if (pType.equals(ESRIGRID)) {
                readArcGrid(mapFile);
            } else if (pType.equals(GEOTIFF)) {
                readGeotiff(mapFile);
            } else if (pType.equals(JPG)) {
                readWorldImage(mapFile);
            } else if (pType.equals(PNG)) {
                readWorldImage(mapFile);
            } else if (pType.equals(AIG) || pType.endsWith("w001001x.adf")) {
                readAig(mapFile);
            } else if (pType.equals(GRASS)) {
                readGrass(mapFile);
            } else {
                throw new ModelsIllegalargumentException("Data type not supported: " + pType, this.getClass().getSimpleName(), pm);
            }
        } finally {
            pm.done();
        }

    }

    private void readGrass( File mapFile ) throws Exception {
        JGrassMapEnvironment mapEnvironment = new JGrassMapEnvironment(new File(file));
        CoordinateReferenceSystem crs = mapEnvironment.getCoordinateReferenceSystem();
        JGrassRegion readRegion = mapEnvironment.getFileRegion();
        double n = readRegion.getNorth();
        double s = readRegion.getSouth();
        double w = readRegion.getWest();
        double e = readRegion.getEast();

        Envelope env = readRegion.getEnvelope();
        originalEnvelope = new GeneralEnvelope(new ReferencedEnvelope(env.getMinX(), env.getMaxX(), env.getMinY(), env.getMaxY(),
                crs));

        // if bounds supplied, use them as region
        if (pBounds != null) {
            // n, s, w, e
            n = pBounds[0];
            s = pBounds[1];
            w = pBounds[2];
            e = pBounds[3];
        }
        if (pRes != null) {
            readRegion = new JGrassRegion(w, e, s, n, pRes[0], pRes[1]);
        }
        if (pRowcol != null) {
            readRegion = new JGrassRegion(w, e, s, n, pRowcol[0], pRowcol[1]);
        }

        if (!doEnvelope) {
            int r = readRegion.getRows();
            int c = readRegion.getCols();
            if (!JGTConstants.doesOverFlow(r, c) && !doLegacyGrassAutoBoxed) {
                if (generalParameter == null) {
                    generalParameter = createGridGeometryGeneralParameter(readRegion.getCols(), readRegion.getRows(),
                            readRegion.getNorth(), readRegion.getSouth(), readRegion.getEast(), readRegion.getWest(), crs);
                }
                GrassCoverageFormat format = new GrassCoverageFormatFactory().createFormat();
                GrassCoverageReader reader = format.getReader(mapEnvironment.getCELL());
                outRaster = (GridCoverage2D) reader.read(generalParameter);
            } else {
                OmsGrassLegacyReader reader = new OmsGrassLegacyReader();
                reader.file = file;
                reader.pm = pm;
                reader.inWindow = GrassLegacyUtilities.jgrassRegion2legacyWindow(readRegion);
                reader.readCoverage();
                outRaster = reader.outGC;
            }
            checkNovalues();
        }
    }

    private void readAig( File mapFile ) throws IllegalArgumentException, IOException {
        final ImageLayout l = new ImageLayout();
        l.setTileGridXOffset(0).setTileGridYOffset(0).setTileHeight(512).setTileWidth(512);

        Hints hints = new Hints();
        hints.add(new RenderingHints(JAI.KEY_IMAGE_LAYOUT, l));

        final URL url = mapFile.toURI().toURL();
        final Object source = url;
        final BaseGDALGridCoverage2DReader reader = new AIGReader(source, hints);
        originalEnvelope = reader.getOriginalEnvelope();
        if (!doEnvelope) {
            GridCoverage2D coverage = (GridCoverage2D) reader.read(generalParameter);
            outRaster = coverage.view(ViewType.GEOPHYSICS);

            resample();
            checkNovalues();
        }
    }

    private void readGeotiff( File mapFile ) throws IOException {
        String nameWithoutExtention = FileUtilities.getNameWithoutExtention(mapFile);
        File twfFile = new File(mapFile.getParentFile(), nameWithoutExtention + ".tfw");
        AbstractGridCoverage2DReader geoTiffReader;
        if (twfFile.exists()) {
            geoTiffReader = new WorldImageReader(mapFile);
        } else {
            geoTiffReader = new GeoTiffReader(mapFile);
        }
        originalEnvelope = geoTiffReader.getOriginalEnvelope();
        if (!doEnvelope) {
            GridCoverage2D coverage = geoTiffReader.read(generalParameter);
            outRaster = coverage.view(ViewType.GEOPHYSICS);

            resample();
            checkNovalues();
        }
    }

    private void readArcGrid( File mapFile ) throws IllegalArgumentException, IOException {
        ArcGridReader arcGridReader = new ArcGridReader(mapFile);
        originalEnvelope = arcGridReader.getOriginalEnvelope();
        if (!doEnvelope) {
            GridCoverage2D coverage = arcGridReader.read(generalParameter);
            outRaster = coverage.view(ViewType.GEOPHYSICS);

            resample();
            checkNovalues();
        }
    }

    private void readWorldImage( File mapFile ) throws IllegalArgumentException, IOException {
        WorldImageReader worldImageReader = new WorldImageReader(mapFile);
        originalEnvelope = worldImageReader.getOriginalEnvelope();
        if (!doEnvelope) {
            GridCoverage2D coverage = worldImageReader.read(generalParameter);
            outRaster = coverage.view(ViewType.RENDERED);

            resample();
            // checkNovalues();
        }
    }

    private void resample() {
        if (!hasBoundsRequest() && !hasResolutionRequest() && !hasRowColsRequest()) {
            // no resample required
            return;
        }

        HashMap<String, Double> envelopeParams = getRegionParamsFromGridCoverage(outRaster);
        double west = envelopeParams.get(WEST);
        double south = envelopeParams.get(SOUTH);
        double east = envelopeParams.get(EAST);
        double north = envelopeParams.get(NORTH);
        double xres = envelopeParams.get(XRES);
        double yres = envelopeParams.get(YRES);

        if (pBounds == null) {
            pBounds = new double[]{north, south, west, east};
        }
        if (pRes == null) {
            pRes = new double[]{xres, yres};
        }

        double n = pBounds[0];
        double s = pBounds[1];
        double w = pBounds[2];
        double e = pBounds[3];
        if (pRes != null || pRowcol != null) {
            int newRows = 0;
            int newCols = 0;
            if (pRowcol != null) {
                newRows = (int) pRowcol[0];
                newCols = (int) pRowcol[1];
                if (pRes == null) {
                    pRes = new double[2];
                }
                pRes[0] = (e - w) / (double) newCols;
                pRes[1] = (n - s) / (double) newRows;
            } else if (pRes != null) {
                pRowcol = new int[2];
                newRows = (int) Math.round((n - s) / pRes[1]);
                newCols = (int) Math.round((e - w) / pRes[0]);
            }
            pRowcol[0] = newRows;
            pRowcol[1] = newCols;
        }

        HashMap<String, Double> newParams = makeRegionParamsMap(n, s, w, e, pRes[0], pRes[1], (int) pRowcol[1], (int) pRowcol[0]);
        CoordinateReferenceSystem crs = outRaster.getCoordinateReferenceSystem();
        GridGeometry2D gg = gridGeometryFromRegionParams(newParams, crs);
        outRaster = (GridCoverage2D) Operations.DEFAULT.resample(outRaster, crs, gg, null);
    }

    private void checkNovalues() {
        // TODO make this nice, this can't be the way
        if (fileNovalue == null || geodataNovalue == null) {
            return;
        }
        if (isNovalue(internalFileNovalue) && isNovalue(internalGeodataNovalue)) {
            return;
        }
        if (!NumericsUtilities.dEq(internalFileNovalue, internalGeodataNovalue)) {
            HashMap<String, Double> params = getRegionParamsFromGridCoverage(outRaster);
            int height = params.get(ROWS).intValue();
            int width = params.get(COLS).intValue();
            WritableRaster tmpWR = createDoubleWritableRaster(width, height, null, null, null);
            WritableRandomIter tmpIter = RandomIterFactory.createWritable(tmpWR, null);
            RenderedImage readRI = outRaster.getRenderedImage();
            RandomIter readIter = RandomIterFactory.create(readRI, null);
            int minX = readRI.getMinX();
            int minY = readRI.getMinY();
            for( int r = 0; r < height; r++ ) {
                for( int c = 0; c < width; c++ ) {
                    double value = readIter.getSampleDouble(c + minX, r + minY, 0);
                    if (isNovalue(value) || value == internalFileNovalue || value == -Float.MAX_VALUE || value == Float.MAX_VALUE) {
                        tmpIter.setSample(c, r, 0, internalGeodataNovalue);
                    } else {
                        tmpIter.setSample(c, r, 0, value);
                    }
                }
            }
            readIter.done();
            tmpIter.done();
            outRaster = buildCoverage(new File(file).getName(), tmpWR, params, outRaster.getCoordinateReferenceSystem());
        }
    }

    private boolean hasBoundsRequest() {
        return pNorth != null && pSouth != null && pWest != null && pEast != null;
    }

    private boolean hasRowColsRequest() {
        return pRows != null && pCols != null;
    }

    private boolean hasResolutionRequest() {
        return pXres != null && pYres != null;
    }

    /**
     * Utility method to quickly read a grid in default mode.
     * 
     * @param path the path to the file.
     * @return the read coverage.
     * @throws Exception
     */
    public static GridCoverage2D readRaster( String path ) throws Exception {
        OmsRasterReader reader = new OmsRasterReader();
        reader.file = path;
        reader.process();
        GridCoverage2D geodata = reader.outRaster;
        return geodata;
    }
}
