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
package org.jgrasstools.gears.io.geopaparazzi.forms;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.jgrasstools.gears.io.geopaparazzi.forms.items.ItemBoolean;
import org.jgrasstools.gears.io.geopaparazzi.forms.items.ItemCombo;
import org.jgrasstools.gears.io.geopaparazzi.forms.items.ItemLabel;
import org.jgrasstools.gears.io.geopaparazzi.forms.items.ItemText;
import org.jgrasstools.gears.io.geopaparazzi.forms.items.ItemTextArea;
import org.jgrasstools.gears.utils.files.FileUtilities;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

/**
 * Form utilities
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class Utilities {
    
    public static final String ATTR_SECTIONNAME = "sectionname";
    public static final String ATTR_SECTIONOBJECTSTR = "sectionobjectstr";
    public static final String ATTR_FORMS = "forms";
    public static final String ATTR_FORMNAME = "formname";

    public static final String TAG_LONGNAME = "longname";
    public static final String TAG_SHORTNAME = "shortname";
    public static final String TAG_FORMS = "forms";
    public static final String TAG_FORMITEMS = "formitems";
    public static final String TAG_KEY = "key";
    public static final String TAG_VALUE = "value";
    public static final String TAG_VALUES = "values";
    public static final String TAG_ITEMS = "items";
    public static final String TAG_ITEM = "item";
    public static final String TAG_TYPE = "type";
    public static final String TAG_READONLY = "readonly";
    public static final String TAG_SIZE = "size";
    public static final String TAG_URL = "url";

    public static void properties2Mainframe( MainFrame mainFrame, File templateFile ) throws Exception {
        List<String> templateLinesList = FileUtilities.readFileToLinesList(templateFile);
        String name = FileUtilities.getNameWithoutExtention(templateFile);

        Section currentSection = new Section(name);
        mainFrame.addSection(currentSection);
        Form currentForm = null;
        for( int i = 0; i < templateLinesList.size(); i++ ) {
            String line = templateLinesList.get(i).trim();
            if (line.length() == 0) {
                continue;
            }
            if (line.startsWith("#")) {
                String title = line.substring(1).trim();
                currentForm = new Form(title);
                currentSection.addForms(currentForm);
                continue;
            }
            String[] split = line.split("\\|");
            String type = split[0].trim();

            if (type.equals("text")) {
                String field = split[1].trim();
                String mandatory = split[2].trim();
                String value = "";
                if (split.length == 4) {
                    value = split[3].trim();
                }

                ItemText item = new ItemText(field, value, Boolean.parseBoolean(mandatory), false);
                currentForm.addItem(item);
            } else if (type.startsWith("textarea")) {
                String field = split[1].trim();
                String mandatory = split[2].trim();
                String value = "";
                if (split.length == 4) {
                    value = split[3].trim();
                }

                ItemTextArea item = new ItemTextArea(field, value, Boolean.parseBoolean(mandatory), false);
                currentForm.addItem(item);
            } else if (type.startsWith("combo")) {
                String field = split[1].trim();
                String mandatory = split[2].trim();
                String value = "";
                if (split.length == 4) {
                    value = split[3].trim();
                }
                String comboItems = type.replaceFirst("combo:", "");
                String[] itemsSplit = comboItems.split(",");
                for( int j = 0; j < itemsSplit.length; j++ ) {
                    itemsSplit[j] = itemsSplit[j].trim();
                }

                ItemCombo combo = new ItemCombo(field, itemsSplit, value, Boolean.parseBoolean(mandatory));
                currentForm.addItem(combo);
            } else if (type.equals("checkbox")) {
                String field = split[1].trim();
                String mandatory = split[2].trim();
                String value = "false";
                if (split.length == 4) {
                    value = split[3].trim();
                }

                ItemBoolean checkbox = new ItemBoolean(field, value, Boolean.parseBoolean(mandatory));
                currentForm.addItem(checkbox);
            } else if (type.equals("label")) {
                String label = "";
                if (split.length > 1)
                    label = split[1].trim();

                ItemLabel labelItem = new ItemLabel(label, 20, false);
                currentForm.addItem(labelItem);
            }
        }
    }
    
    public static List<String> getFormNames4Section( JSONObject section ) throws JSONException {
        List<String> names = new ArrayList<String>();
        JSONArray jsonArray = section.getJSONArray(ATTR_FORMS);
        if (jsonArray != null && jsonArray.length() > 0) {
            for( int i = 0; i < jsonArray.length(); i++ ) {
                JSONObject jsonObject = jsonArray.getJSONObject(i);
                if (jsonObject.has(ATTR_FORMNAME)) {
                    String formName = jsonObject.getString(ATTR_FORMNAME);
                    names.add(formName);
                }
            }
        }
        return names;
    }

    public static JSONObject getForm4Name( String formName, JSONObject section ) throws JSONException {
        JSONArray jsonArray = section.getJSONArray(ATTR_FORMS);
        if (jsonArray != null && jsonArray.length() > 0) {
            for( int i = 0; i < jsonArray.length(); i++ ) {
                JSONObject jsonObject = jsonArray.getJSONObject(i);
                if (jsonObject.has(ATTR_FORMNAME)) {
                    String tmpFormName = jsonObject.getString(ATTR_FORMNAME);
                    if (tmpFormName.equals(formName)) {
                        return jsonObject;
                    }
                }
            }
        }
        return null;
    }
    
    /**
     * Utility method to get the formitems of a form object.
     * 
     * <p>Note that the entering json object has to be one 
     * object of the main array, not THE main array itself, 
     * i.e. a choice was already done.
     * 
     * @param jsonObj the single object.
     * @return the array of items of the contained form or <code>null</code> if 
     *          no form is contained.
     * @throws JSONException
     */
    public static JSONArray getFormItems( JSONObject formObj ) throws JSONException {
        if (formObj.has(TAG_FORMITEMS)) {
            JSONArray formItemsArray = formObj.getJSONArray(TAG_FORMITEMS);
            return formItemsArray;
        }
        return null;
    }

}
