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
package org.jgrasstools.hortonmachine.models.hm;

import java.util.HashMap;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.AdigeModified;
import org.jgrasstools.hortonmachine.utils.HMTestCase;

public class TestNewAge_runoff_2 extends HMTestCase {
	public static final double doubleNovalue = Double.NaN;
	
     public void testAb() throws Exception {
    	 
		String startDate = "1996-01-01 00:00";
		String endDate = "1996-01-01 00:00";
		int timeStepMinutes = 60;
	//	String fId = "ID";
		
		OmsShapefileFeatureReader networkreader = new OmsShapefileFeatureReader();
		networkreader.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/pfaf_tca400_mode2_6.shp";
		networkreader.readFeatureCollection();
		SimpleFeatureCollection networkFC = networkreader.geodata;
		
		
		OmsShapefileFeatureReader subbasinreader = new OmsShapefileFeatureReader();
		subbasinreader.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/subbasin_tca4000_mode2_6.shp";
		subbasinreader.readFeatureCollection();
		SimpleFeatureCollection subbasinreaderFC = subbasinreader.geodata;
		
		OmsShapefileFeatureReader hydrometerreader = new OmsShapefileFeatureReader();
		hydrometerreader.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/hydro4.shp";
		hydrometerreader.readFeatureCollection();
		SimpleFeatureCollection hydrometerreaderFC = hydrometerreader.geodata;
		

		OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
		reader_rain.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/output/P_DK_interpolated.csv";
		reader_rain.idfield = "ID";
		reader_rain.tStart = startDate;
		reader_rain.tEnd = endDate;
		reader_rain.fileNovalue = "-999.0";
		reader_rain.tTimestep = timeStepMinutes;
	//	reader_rain.initProcess();

    	OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
		reader_hydro.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/Hydrometers1994_2006_22.csv";
		reader_hydro.idfield = "ID";
		reader_hydro.tStart = startDate;
		reader_hydro.tEnd = endDate;
		reader_hydro.fileNovalue = "-9999";
		reader_hydro.tTimestep = timeStepMinutes;
	//	reader_hydro.initProcess();

		OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
		reader_et.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/krigingg_interpolated_etp2.csv";
		reader_et.idfield = "ID";
		reader_et.tStart = startDate;
		reader_et.tEnd = endDate;
		reader_et.fileNovalue = "-9999";
		reader_et.tTimestep = 60;
		//reader_et.initProcess();
		
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/output/discharge_test.csv";
		writer.tStart = reader_rain.tStart;
		writer.tTimestep = reader_rain.tTimestep;
		writer.fileNovalue = "-9999";
		writer.tTimestep = timeStepMinutes;
		
		
    AdigeModified adige = new AdigeModified();
        adige.pm = pm;
    	adige.pCmax = 106.55212488441822;
		adige.pB = 0.971346077530604;
		adige.pAlpha = 0.8084184792349338;
		adige.pRs = 0.024261525315367;
		adige.pRq = 0.0040606747636673;
		adige.pQ0 = 0.1;

		
     	adige.pmode = 1;
    	adige.pDorouting = false;
     //	adige.pDoGeom = false;
     	//adige.fIdDistanze = "dtoou";
    	adige.fIdnetnum = "netnum";
     	adige.inHillslope = subbasinreaderFC;
     	adige.fNetnum = "netnum";
     	adige.fBaricenter = "avgZ";
     	// adige.fVegetation = "uso_reclas";
     	adige.inHydrometers = hydrometerreaderFC;
     	// adige.inDams = damsFC;
     	// adige.inTributary = tributaryFC;
     	// adige.inOfftakes = offtakesFC;
     	// adige.inVegetation = vegetationData;
     	adige.pPfafids = "1";
     	adige.fMonpointid = "ID";
     	adige.inNetwork = networkFC;
     	adige.fPfaff = "pfafstette";
     	adige.fNetelevstart = "elevfirstp";
     	adige.fNetelevend = "elevlastpo";
     	adige.pRainintensity = -1;
     	adige.pRainduration = -1;
     	adige.doLog = false;
        adige.pNetNumCali = 9;
	
		adige.pDoReadMeas = false;
		adige.pPathtoMeas = "/Users/administrator/Desktop/test_dischareg";
		adige.pDimMeas = 25;
		
    	//adige.tTimestep = timeStepMinutes;
		adige.tStart = startDate;
		adige.tEnd = endDate;
		reader_rain.initProcess();
		
		int indmod = 0;
     while(reader_rain.doProcess)  {
			reader_rain.nextRecord();
			HashMap<Integer, double[]> rainvalues = reader_rain.outData;
			reader_et.nextRecord();
			HashMap<Integer, double[]> etvalues = reader_et.outData;
         //  reader_melting.nextRecord();
         //   HashMap<Integer, double[]>  meltingvalues = reader_melting.outData;
			adige.inRain = rainvalues;
			adige.inEtp = etvalues;
         // adige.inMelting = meltingvalues;
         //  
			reader_hydro.nextRecord();
          adige.inHydrometerdata = reader_hydro.outData;
	      
          adige.process();
            

 
			HashMap<Integer, double[]> outdischarge = adige.outDischarge;
	        writer.inData = outdischarge;
            writer.writeNextLine();
            
        	indmod++;
            
        
        }
   
     
     System.out.println("ciao");
        reader_rain.close();
        reader_et.close();
        reader_hydro.close();
        writer.close();
        
      
    
}


}
