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

/**
 * A boolean item.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class ItemBoolean implements Item {

    private String description;
    private String defaultValue;
    private boolean isMandatory;

    /**
     * @param description
     * @param defaultValue a default value: <b>false</b> or <b>true</b>.
     * @param isMandatory
     */
    public ItemBoolean( String description, String defaultValue, boolean isMandatory ) {
        this.isMandatory = isMandatory;
        if (defaultValue == null) {
            defaultValue = "false";
        }
        this.description = description;
        this.defaultValue = defaultValue;
        checkBoolean(this.defaultValue);
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("        {\n");
        sb.append("             \"key\": \"").append(description).append("\",\n");
        sb.append("             \"value\": \"").append(defaultValue).append("\",\n");
        sb.append("             \"type\": \"").append("boolean").append("\",\n");
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
        checkBoolean(value);
        defaultValue = value;
    }

    private void checkBoolean( String value ) {
        if (!value.equals("false") && !value.equals("true")) {
            throw new IllegalArgumentException("Value has to be false or true.");
        }
    }

    @Override
    public String getValue() {
        return defaultValue;
    }
}
