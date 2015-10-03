/*
 * JGrass - Free Open Source Java GIS http://www.jgrass.org 
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Library General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Library General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU Library General Public License
 * along with this library; if not, write to the Free Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
package org.jgrasstools.hortonmachine.modules.networktools.epanet.core.types;

import org.joda.time.DateTime;

/**
 * A class representing pipes data for a given time instant.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class Pipe {
    /**
     * Current time instant.
     */
    public DateTime time;

    public String id;

    public float[] flow = new float[2];

    public float[] velocity = new float[2];

    public float headloss;

    public float status;

    @SuppressWarnings("nls")
    @Override
    public String toString() {
        return id + "\t" + flow[0] + "\t" + flow[1] + "\t" + velocity[0] + "\t" + velocity[1] + "\t" + headloss + "\t" + status;
    }

}
