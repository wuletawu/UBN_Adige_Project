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
package org.jgrasstools.gears.io.geopaparazzi.forms.items;

import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * A date item.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class ItemDate implements Item {

    private String description;
    private boolean isMandatory;
    private String defaultValueStr;

    private SimpleDateFormat dateFormatter = new SimpleDateFormat("yyyy-MM-dd");

    public ItemDate( String description, Date defaultValue, boolean isMandatory ) {
        if (defaultValue == null) {
            defaultValueStr = "";
        } else {
            this.defaultValueStr = dateFormatter.format(defaultValue);
        }
        this.description = description;
        this.isMandatory = isMandatory;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("        {\n");
        sb.append("             \"key\": \"").append(description).append("\",\n");
        sb.append("             \"value\": \"").append(defaultValueStr).append("\",\n");
        sb.append("             \"type\": \"").append("date").append("\",\n");
        sb.append("             \"mandatory\": \"").append(isMandatory ? "yes" : "no").append("\"\n");
        sb.append("        }\n");
        return sb.toString();
    }

    @Override
    public String getKey() {
        return description;
    }

    @Override
    public void setValue( String value ) {
        defaultValueStr = value;
    }

    @Override
    public String getValue() {
        return defaultValueStr;
    }
}
