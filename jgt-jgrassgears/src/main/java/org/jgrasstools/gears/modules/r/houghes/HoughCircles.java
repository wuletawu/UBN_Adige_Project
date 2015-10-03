package org.jgrasstools.gears.modules.r.houghes;
/** houghCircles_.java:
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 @author Hemerson Pistori (pistori@ec.ucdb.br) and Eduardo Rocha Costa
 @created 18 de Mar�o de 2004
 
 The Hough Transform implementation was based on 
 Mark A. Schulze applet (http://www.markschulze.net/)
 
*/

//package sigus.templateMatching;
//import sigus.*;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.imageio.ImageIO;
import javax.media.jai.iterator.RandomIter;
import javax.media.jai.iterator.RandomIterFactory;

import com.vividsolutions.jts.geom.Coordinate;

/**
 * A Hough Transform implementation.
 * 
 * @author Hemerson Pistori (pistori@ec.ucdb.br) and Eduardo Rocha Costa
 * @author Mark A. Schulze applet (http://www.markschulze.net/)
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class HoughCircles {

    public int radiusMin; // Find circles with radius grater or equal radiusMin
    public int radiusMax; // Find circles with radius less or equal radiusMax
    public int radiusInc; // Increment used to go from radiusMin to radiusMax
    public int maxCircles; // Numbers of circles to be found
    public int threshold = -1; // An alternative to maxCircles. All circles with
    // a value in the hough space greater then threshold are marked. Higher thresholds
    // results in fewer circles.
    byte imageValues[]; // Raw image (returned by ip.getPixels())
    public int width; // Hough Space width (depends on image width)
    public int height; // Hough Space heigh (depends on image height)
    public int depth; // Hough Space depth (depends on radius interval)
    public int offset; // Image Width
    public int offx; // ROI x offset
    public int offy; // ROI y offset
    private int vectorMaxSize = 500;
    boolean useThreshold = false;
    int lut[][][]; // LookUp Table for rsin e rcos values
    private BufferedImage raster;

    public static void main( String[] args ) throws Exception {

        int[] i = {10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
        for( int index : i ) {
            String inImage = "/home/hydrologis/data/rilievo_tls/slice_" + index + ".0rgb.png";
            String outImage = "/home/hydrologis/data/rilievo_tls/slice_" + index + ".0bw_out.png";

            BufferedImage src = ImageIO.read(new File(inImage));
            HoughCircles h = new HoughCircles(src, 10, 40, 1, 300);
            h.run();
            ImageIO.write(src, "png", new File(outImage));
        }

    }

    public HoughCircles( BufferedImage raster, int radiusMin, int radiusMax, int radiusIncrement, int circleCount ) {
        this.raster = raster;
        this.radiusMin = radiusMin;
        this.radiusMax = radiusMax;
        maxCircles = circleCount;
        radiusInc = radiusIncrement;
        depth = ((radiusMax - radiusMin) / radiusInc) + 1;

        // TODO support threshold
        // threshold = (int) gd.getNextNumber();
        // if (maxCircles > 0) {
        // useThreshold = false;
        // threshold = -1;
        // } else {
        // useThreshold = true;
        // if (threshold < 0) {
        // IJ.showMessage("Threshold must be greater than 0");
        // return (false);
        // }
        // }
    }

    public void run() {
        offx = 0;
        offy = 0;
        width = raster.getWidth();
        height = raster.getHeight();
        offset = width;

        imageValues = new byte[width * height];
        int count = 0;
        RandomIter renderedImageIterator = RandomIterFactory.create(raster, null);
        for( int r = 0; r < height; r++ ) {
            for( int c = 0; c < width; c++ ) {
                int sample = renderedImageIterator.getSample(c, r, 0);
                imageValues[count++] = (byte) sample;
            }
        }
        renderedImageIterator.done();

        double[][][] houghValues = houghTransform();

        // Create image View for Hough Transform.
        // ImageProcessor newip = new ByteProcessor(width, height);
        // byte[] newpixels = (byte[]) newip.getPixels();
        // createHoughPixels(houghValues, newpixels);
        //
        // // Create image View for Marked Circles.
        // ImageProcessor circlesip = new ByteProcessor(width, height);
        // byte[] circlespixels = (byte[]) circlesip.getPixels();

        // Mark the center of the found circles in a new image
        // if (useThreshold)
        // getCenterPointsByThreshold(threshold);
        // else
        Coordinate[] centerPoints = getCenterPoints(houghValues, maxCircles);
        // drawCircles(houghValues, circlespixels);

        Graphics2D g2d = (Graphics2D) raster.getGraphics();
        g2d.setColor(Color.red);
        g2d.setStroke(new BasicStroke(2));
        for( Coordinate point : centerPoints ) {
            int size = (int) point.z * 2;
            g2d.drawOval((int) point.x - size / 2, (int) point.y + -size / 2, size, size);
        }
    }

    /** The parametric equation for a circle centered at (a,b) with
        radius r is:

    a = x - r*cos(theta)
    b = y - r*sin(theta)

    In order to speed calculations, we first construct a lookup
    table (lut) containing the rcos(theta) and rsin(theta) values, for
    theta varying from 0 to 2*PI with increments equal to
    1/8*r. As of now, a fixed increment is being used for all
    different radius (1/8*radiusMin). This should be corrected in
    the future.

    Return value = Number of angles for each radius
       
    */
    private int buildLookUpTable() {

        int i = 0;
        int incDen = Math.round(8F * radiusMin); // increment denominator

        lut = new int[2][incDen][depth];

        for( int radius = radiusMin; radius <= radiusMax; radius = radius + radiusInc ) {
            i = 0;
            for( int incNun = 0; incNun < incDen; incNun++ ) {
                double angle = (2 * Math.PI * (double) incNun) / (double) incDen;
                int indexR = (radius - radiusMin) / radiusInc;
                int rcos = (int) Math.round((double) radius * Math.cos(angle));
                int rsin = (int) Math.round((double) radius * Math.sin(angle));
                if ((i == 0) | (rcos != lut[0][i][indexR]) & (rsin != lut[1][i][indexR])) {
                    lut[0][i][indexR] = rcos;
                    lut[1][i][indexR] = rsin;
                    i++;
                }
            }
        }

        return i;
    }

    private double[][][] houghTransform() {
        int lutSize = buildLookUpTable();
        double[][][] houghValues = new double[width][height][depth];
        int k = width - 1;
        int l = height - 1;

        for( int y = 1; y < l; y++ ) {
            for( int x = 1; x < k; x++ ) {
                for( int radius = radiusMin; radius <= radiusMax; radius = radius + radiusInc ) {
                    if (imageValues[(x + offx) + (y + offy) * offset] != 0) {// Edge pixel found
                        int indexR = (radius - radiusMin) / radiusInc;
                        for( int i = 0; i < lutSize; i++ ) {
                            int a = x + lut[1][i][indexR];
                            int b = y + lut[0][i][indexR];
                            if ((b >= 0) & (b < height) & (a >= 0) & (a < width)) {
                                houghValues[a][b][indexR] += 1;
                            }
                        }
                    }
                }
            }
        }

        return houghValues;

    }

    // Convert Values in Hough Space to an 8-Bit Image Space.
    private void createHoughPixels( double[][][] houghValues, byte houghPixels[] ) {
        double d = -1D;
        for( int j = 0; j < height; j++ ) {
            for( int k = 0; k < width; k++ ) {
                if (houghValues[k][j][0] > d) {
                    d = houghValues[k][j][0];
                }
            }
        }
        for( int l = 0; l < height; l++ ) {
            for( int i = 0; i < width; i++ ) {
                houghPixels[i + l * width] = (byte) Math.round((houghValues[i][l][0] * 255D) / d);
            }
        }
    }

    // Draw the circles found in the original image.
    public void drawCircles( double[][][] houghValues, byte[] circlespixels ) {

        // Copy original input pixels into output
        // circle location display image and
        // combine with saturation at 100
        int roiaddr = 0;
        for( int y = offy; y < offy + height; y++ ) {
            for( int x = offx; x < offx + width; x++ ) {
                // Copy;
                circlespixels[roiaddr] = imageValues[x + offset * y];
                // Saturate
                if (circlespixels[roiaddr] != 0)
                    circlespixels[roiaddr] = 100;
                else
                    circlespixels[roiaddr] = 0;
                roiaddr++;
            }
        }
        // Copy original image to the circlespixels image.
        // Changing pixels values to 100, so that the marked
        // circles appears more clear. Must be improved in
        // the future to show the resuls in a colored image.
        // for(int i = 0; i < width*height ;++i ) {
        // if(imageValues[i] != 0 )
        // if(circlespixels[i] != 0 )
        // circlespixels[i] = 100;
        // else
        // circlespixels[i] = 0;
        // }
        // if (centerPoints == null) {
        // if (useThreshold)
        // getCenterPointsByThreshold(threshold);
        // else
        Coordinate[] centerPoints = getCenterPoints(houghValues, maxCircles);
        // }
        byte cor = -1;
        // Redefine these so refer to ROI coordinates exclusively
        int offset = width;
        int offx = 0;
        int offy = 0;

        for( int l = 0; l < maxCircles; l++ ) {
            int i = (int) centerPoints[l].x;
            int j = (int) centerPoints[l].y;
            // Draw a gray cross marking the center of each circle.
            for( int k = -10; k <= 10; ++k ) {
                int p = (j + k + offy) * offset + (i + offx);
                if (!outOfBounds(j + k + offy, i + offx))
                    circlespixels[(j + k + offy) * offset + (i + offx)] = cor;
                if (!outOfBounds(j + offy, i + k + offx))
                    circlespixels[(j + offy) * offset + (i + k + offx)] = cor;
            }
            for( int k = -2; k <= 2; ++k ) {
                if (!outOfBounds(j - 2 + offy, i + k + offx))
                    circlespixels[(j - 2 + offy) * offset + (i + k + offx)] = cor;
                if (!outOfBounds(j + 2 + offy, i + k + offx))
                    circlespixels[(j + 2 + offy) * offset + (i + k + offx)] = cor;
                if (!outOfBounds(j + k + offy, i - 2 + offx))
                    circlespixels[(j + k + offy) * offset + (i - 2 + offx)] = cor;
                if (!outOfBounds(j + k + offy, i + 2 + offx))
                    circlespixels[(j + k + offy) * offset + (i + 2 + offx)] = cor;
            }
        }
    }

    private boolean outOfBounds( int y, int x ) {
        if (x >= width)
            return (true);
        if (x <= 0)
            return (true);
        if (y >= height)
            return (true);
        if (y <= 0)
            return (true);
        return (false);
    }

    /** Search for a fixed number of circles.

    @param maxCircles The number of circles that should be found.  
     * @param houghValues 
    */
    private Coordinate[] getCenterPoints( double[][][] houghValues, int maxCircles ) {

        Coordinate[] centerPoints = new Coordinate[maxCircles];
        int xMax = 0;
        int yMax = 0;
        int rMax = 0;

        for( int c = 0; c < maxCircles; c++ ) {
            double counterMax = -1;
            for( int radius = radiusMin; radius <= radiusMax; radius = radius + radiusInc ) {

                int indexR = (radius - radiusMin) / radiusInc;
                for( int y = 0; y < height; y++ ) {
                    for( int x = 0; x < width; x++ ) {
                        if (houghValues[x][y][indexR] > counterMax) {
                            counterMax = houghValues[x][y][indexR];
                            xMax = x;
                            yMax = y;
                            rMax = radius;
                        }
                    }

                }
            }

            centerPoints[c] = new Coordinate(xMax, yMax, rMax);

            clearNeighbours(houghValues, xMax, yMax, rMax);
        }
        return centerPoints;
    }

    /** Search circles having values in the hough space higher than a threshold

    @param threshold The threshold used to select the higher point of Hough Space
    */
    // private void getCenterPointsByThreshold( int threshold ) {
    //
    // centerPoint = new Point[vectorMaxSize];
    // int xMax = 0;
    // int yMax = 0;
    // int countCircles = 0;
    //
    // for( int radius = radiusMin; radius <= radiusMax; radius = radius + radiusInc ) {
    // int indexR = (radius - radiusMin) / radiusInc;
    // for( int y = 0; y < height; y++ ) {
    // for( int x = 0; x < width; x++ ) {
    //
    // if (houghValues[x][y][indexR] > threshold) {
    //
    // if (countCircles < vectorMaxSize) {
    //
    // centerPoint[countCircles] = new Point(x, y);
    //
    // clearNeighbours(xMax, yMax, radius);
    //
    // ++countCircles;
    // } else
    // break;
    // }
    // }
    // }
    // }
    //
    // maxCircles = countCircles;
    // }

    /** Clear, from the Hough Space, all the counter that are near (radius/2) a previously found circle C.
        
    @param x The x coordinate of the circle C found.
    @param x The y coordinate of the circle C found.
    @param x The radius of the circle C found.
    */
    private void clearNeighbours( double[][][] houghValues, int x, int y, int radius ) {

        // The following code just clean the points around the center of the circle found.

        double halfRadius = radius / 2.0F;
        double halfSquared = halfRadius * halfRadius;

        int y1 = (int) Math.floor((double) y - halfRadius);
        int y2 = (int) Math.ceil((double) y + halfRadius) + 1;
        int x1 = (int) Math.floor((double) x - halfRadius);
        int x2 = (int) Math.ceil((double) x + halfRadius) + 1;

        if (y1 < 0)
            y1 = 0;
        if (y2 > height)
            y2 = height;
        if (x1 < 0)
            x1 = 0;
        if (x2 > width)
            x2 = width;

        for( int r = radiusMin; r <= radiusMax; r = r + radiusInc ) {
            int indexR = (r - radiusMin) / radiusInc;
            for( int i = y1; i < y2; i++ ) {
                for( int j = x1; j < x2; j++ ) {
                    if (Math.pow(j - x, 2D) + Math.pow(i - y, 2D) < halfSquared) {
                        houghValues[j][i][indexR] = 0.0D;
                    }
                }
            }
        }

    }

}
