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
package org.jgrasstools.gears.io.eicalculator;

import static org.jgrasstools.gears.i18n.GearsMessages.OMSEIENERGYREADER_AUTHORCONTACTS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSEIENERGYREADER_AUTHORNAMES;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSEIENERGYREADER_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSEIENERGYREADER_KEYWORDS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSEIENERGYREADER_LABEL;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSEIENERGYREADER_LICENSE;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSEIENERGYREADER_NAME;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSEIENERGYREADER_STATUS;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSEIENERGYREADER_file_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSEIENERGYREADER_outEnergy_DESCRIPTION;
import static org.jgrasstools.gears.i18n.GearsMessages.OMSEIENERGYREADER_pSeparator_DESCRIPTION;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.Finalize;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;
import oms3.annotations.UI;

import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;

@Description(OMSEIENERGYREADER_DESCRIPTION)
@Author(name = OMSEIENERGYREADER_AUTHORNAMES, contact = OMSEIENERGYREADER_AUTHORCONTACTS)
@Keywords(OMSEIENERGYREADER_KEYWORDS)
@Label(OMSEIENERGYREADER_LABEL)
@Name(OMSEIENERGYREADER_NAME)
@Status(OMSEIENERGYREADER_STATUS)
@License(OMSEIENERGYREADER_LICENSE)
public class OmsEIEnergyReader extends JGTModel {

    @Description(OMSEIENERGYREADER_file_DESCRIPTION)
    @UI(JGTConstants.FILEIN_UI_HINT)
    @In
    public String file = null;

    @Description(OMSEIENERGYREADER_pSeparator_DESCRIPTION)
    @In
    public String pSeparator = ",";

    @Description(OMSEIENERGYREADER_outEnergy_DESCRIPTION)
    @Out
    public List<EIEnergy> outEnergy;

    private BufferedReader csvReader;

    private void ensureOpen() throws IOException {
        if (csvReader == null)
            csvReader = new BufferedReader(new FileReader(file));
    }

    @Finalize
    public void close() throws IOException {
        csvReader.close();
    }

    @Execute
    public void read() throws IOException {
        if (!concatOr(outEnergy == null, doReset)) {
            return;
        }
        ensureOpen();
        outEnergy = new ArrayList<EIEnergy>();
        String line = null;
        while( (line = csvReader.readLine()) != null ) {
            if (line.trim().length() == 0 || line.trim().startsWith("#")) {
                // jump empty lines and lines that start as comment
                continue;
            }
            String[] lineSplit = line.split(pSeparator);
            if (lineSplit.length > 4) {
                throw new IOException("Energy values are defined in 4 columns.");
            }

            EIEnergy eiEnergy = new EIEnergy();
            eiEnergy.basinId = Integer.parseInt(lineSplit[0].trim());
            eiEnergy.energeticBandId = Integer.parseInt(lineSplit[1].trim());
            eiEnergy.virtualMonth = Integer.parseInt(lineSplit[2].trim());
            eiEnergy.energyValue = Double.parseDouble(lineSplit[3].trim());
            outEnergy.add(eiEnergy);
        }

    }

}
