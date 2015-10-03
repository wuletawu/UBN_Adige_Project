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
package org.jgrasstools.gears.io.las.core;

import org.geotools.geometry.jts.ReferencedEnvelope3D;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

/**
 * Las header interface. 
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public interface ILasHeader {

    /**
     * @return the version of the las file as major.minor.
     */
    public abstract String getVersion();

    /**
     * @return the {@link CoordinateReferenceSystem crs} of the file.
     */
    public abstract CoordinateReferenceSystem getCrs();

    /**
     * @return the 3D data envelope of the file.
     */
    public abstract ReferencedEnvelope3D getDataEnvelope();

    /**
     * @return the number of records.
     */
    public long getRecordsCount();

    public abstract boolean hasGpsTime();

    public abstract boolean hasRGB();

    public abstract byte getPointDataFormat();

    public abstract int getGpsTimeType();

    long getOffset();

    short getRecordLength();

    double[] getXYZScale();

    double[] getXYZOffset();
}