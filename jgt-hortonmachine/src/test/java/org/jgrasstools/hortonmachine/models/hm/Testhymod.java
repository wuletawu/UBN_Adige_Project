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
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.AdigeModified_2;
import org.jgrasstools.hortonmachine.utils.HMTestCase;

public class Testhymod extends   HMTestCase {
	public static final double doubleNovalue = Double.NaN;
	
     public void testhymod() throws Exception {
    	 
		String startDate = "2000-01-01 00:00";
		String endDate = "2000-01-10 05:00";
		int timeStepMinutes = 60;
		String fId = "ID";
		
		
		
		OmsShapefileFeatureReader networkreader = new OmsShapefileFeatureReader();
		networkreader.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/pfaf_tca400_mode2_6.shp";
		networkreader.readFeatureCollection();
		networkreader.pm = pm;
		SimpleFeatureCollection networkFC = networkreader.geodata;
		
		
		OmsShapefileFeatureReader subbasinreader = new OmsShapefileFeatureReader();
		subbasinreader.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/subbasin_tca4000_mode2_6.shp";
		subbasinreader.readFeatureCollection();
		subbasinreader.pm = pm;
		SimpleFeatureCollection subbasinreaderFC = subbasinreader.geodata;
		
		OmsShapefileFeatureReader hydrometerreader = new OmsShapefileFeatureReader();
		hydrometerreader.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/hydro4.shp";
		hydrometerreader.readFeatureCollection();
		hydrometerreader.pm = pm;
		SimpleFeatureCollection hydrometerreaderFC = hydrometerreader.geodata;
		

		OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
		reader_rain.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/output/P_DK_interpolated.csv";
		reader_rain.idfield = "ID";
		reader_rain.tStart = startDate;
		reader_rain.tEnd = endDate;
		reader_rain.fileNovalue = "-999.0";
		reader_rain.tTimestep = timeStepMinutes;
		reader_rain.pm = pm;
		
		//reader_rain.initProcess();
		

    	OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
		reader_hydro.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/Hydrometers1995_2013_Novalue.csv";
		reader_hydro.idfield = "ID";
		reader_hydro.tStart = reader_rain.tStart;
		reader_hydro.tEnd = reader_rain.tEnd;
		reader_hydro.fileNovalue = "-9999";
		reader_hydro.tTimestep = reader_rain.tTimestep;
		reader_hydro.pm = pm;
	//	reader_hydro.initProcess();

		OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
		reader_et.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/krigingg_interpolated_etp2.csv";
		reader_et.idfield = "ID";
		reader_et.tStart = reader_rain.tStart;
		reader_et.tEnd = reader_rain.tEnd;
		reader_et.fileNovalue = "-9999";
		reader_et.tTimestep = reader_rain.tTimestep;
		reader_et.pm=pm;
		//reader_et.initProcess();
		
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
	//	TimeSeriesIteratorWriter writer = new TimeSeriesIteratorWriter()
		writer.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/output/discharge_test.csv";
		writer.tStart = reader_rain.tStart;
		writer.tTimestep = reader_rain.tTimestep;
	//	writer.fileNovalue = "-9999";
		writer.tTimestep = reader_rain.tTimestep;
		

		
        AdigeModified_2 adige = new AdigeModified_2();
        adige.pm = pm;
    	adige.pCmax = 106.55212488441822;
		adige.pB = 0.971346077530604;
		adige.pAlpha = 0.8084184792349338;
		adige.pRs = 0.024261525315367;
		adige.pRq = 0.0040606747636673;
		adige.pQ0 = 0.1;
		adige.pMetp = 1.01;
		adige.pMrain = 1.2;

		
     	adige.pmode = 1;
    	adige.pDorouting = false;
     	adige.pDoGeom = true;
     	adige.fIdDistanze = "dtoou";
    	adige.fIdnetnum = "netnum";
     	adige.inHillslope = subbasinreaderFC;
     	adige.inSubbasinDist = subbasinreaderFC;
     	adige.fNetnum = "netnum";
     	adige.fBaricenter = "avgZ";
     	// adige.fVegetation = "uso_reclas";
     	adige.inHydrometers = hydrometerreaderFC;
     	// adige.inDams = damsFC;
     	// adige.inTributary = tributaryFC;
     	// adige.inOfftakes = offtakesFC;
     	// adige.inVegetation = vegetationData;
     	adige.pPfafids = "1,4.5";
     	adige.fMonpointid = "ID";
     	adige.inNetwork = networkFC;
     	adige.fPfaff = "pfafstette";
     	adige.fNetelevstart = "elevfirstp";
     	adige.fNetelevend = "elevlastpo";
     	adige.pRainintensity = -1;
     	adige.pRainduration = -1;
     	adige.doLog = false;
    //    adige.pNetNumCali = 9;
	
		adige.tStart = reader_rain.tStart;
		adige.tEnd = reader_rain.tEnd;
		adige.tTimestep = reader_rain.tTimestep;
		
		     int indmod = 0;
			/*adige.pDoReadMeas = false;
			adige.pPathtoMeas = "/Users/administrator/Desktop/test_dischareg";
			adige.pDimMeas = 25;
			adige.pNetNumCali = 9;*/
		
		reader_rain.initProcess();
     while(reader_rain.doProcess)  {
			reader_rain.nextRecord();
			//System.out.println(reader_rain.tCurrent);
			HashMap<Integer, double[]> rainvalues = reader_rain.outData;
			reader_et.nextRecord();
			HashMap<Integer, double[]> etvalues = reader_et.outData;
         //  reader_melting.nextRecord();
         //   HashMap<Integer, double[]>  meltingvalues = reader_melting.outData;
			adige.inRain = rainvalues;
			adige.inEtp = etvalues;
         // adige.inMelting = meltingvalues;
         //  
		//	reader_hydro.nextRecord();
        //   adige.inHydrometerdata = reader_hydro.outData;
		
			adige.pm = pm;
			adige.process();
         
			
			indmod++;
 
			HashMap<Integer, double[]> outdischarge = adige.outDischarge;
	        writer.inData = outdischarge;
            writer.writeNextLine();
         
            
        
        }
   
   
        reader_rain.close();
        reader_et.close();
       // reader_hydro.close();
        writer.close();
        System.out.println("ciao");    
      
    
}




}
