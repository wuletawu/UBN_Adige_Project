/*package org.jgrasstools.hortonmachine.models.hm;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.libs.monitor.PrintStreamProgressMonitor;
import org.jgrasstools.gears.utils.HMTestCase;
@SuppressWarnings("nls")
public class TestNewAge_runoff_RRS extends HMTestCase {

public TestNewAge_runoff_RRS() throws Exception {
	PrintStreamProgressMonitor pm = new PrintStreamProgressMonitor(System.out, System.out);
	
		String startDate = "1996-01-01 00:00";
		String endDate = "1996-01-10 00:00";
		//int timeStepMinutes = 60;
		
		
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
		reader_rain.fileNovalue = "-9999.0";
		reader_rain.tTimestep = 60;
		reader_rain.initProcess();

		OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
		reader_hydro.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/Hydrometers1994_2006_22.csv";
		reader_hydro.idfield = "ID";
		reader_hydro.tStart = startDate;
		reader_hydro.tEnd = endDate;
		reader_hydro.fileNovalue = "-9999";
		reader_hydro.tTimestep = timeStepMinutes;
		reader_hydro.initProcess();

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
		writer.tStart = startDate;
		//writer.tend = endDate;
		writer.fileNovalue = "-9999";
		writer.tTimestep = 60;
		
		
        RRS adige = new RRS();
        adige.pCmax = 500.0;
     	adige.pB = 2.0;
     	adige.pAlpha = 0.5;
     	adige.pRs = 0.0001;
     	adige.pRq = 0.1;
     	adige.pCff = 10;
     	adige.pCmf = 0.1;
     	adige.pCr = 3;
     	adige.pCs = 4;
        adige.pQ0 = 1.01;
        adige.doDaily = false;
        adige.fBaricenter = "avgZ";
        adige.fMonpointid = "ID";
        adige.fNetelevstart = "elevfirstp";
     	adige.fNetelevend = "elevlastpo";
     	adige.fNetnum = "netnum";
     	adige.fPfaff =  "pfafstette";
     	adige.pathToAirTemperature = "/Users/administrator/Documents/posina/Temp.prj.newage/output/other/Subbasin_LDK_GAUS_temp.csv";
     	adige.pathToDem = "/Users/administrator/Documents/posina/RRSNOW.prj.newage/data/DEM.asc";
     	adige.pathToDischargeOutput = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/output/discharge_test.csv";
     	adige.pathToEtp = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/krigingg_interpolated_etp2.csv";
     	adige.pathToHydrometers = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/Hydrometers1994_2006_22.csv";
     //	adige.pathToHydrometersShape = 
     //	adige.pathToMeltingOutput = 
     	adige.pathToNetworkShape = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/pfaf_tca400_mode2_6.shp";
     	adige.pathToRainfall = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/output/P_DK_interpolated.csv";
     	adige.pathToSkyview = "/Users/administrator/Desktop/RRSNOW_editable.prj.newage/data/skyview.asc";
     	adige.pathToSubBasinsCentroidsShape = "/Users/administrator/Documents/posina/RRSNOW.prj.newage/data/centroid.shp";
     	adige.centroidsID = "netnum";
     	adige.pathToSubbasinShape = subbasinreader.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/subbasin_tca4000_mode2_6.shp";
     	adige.centroidsID = "netnum";
     	//adige.inHillslope = subbasinreaderFC;
     	adige.fNetnum = "netnum";
     	adige.fBaricenter = "avgZ";
     	// adige.fVegetation = "uso_reclas";
     //	adige.inHydrometers = hydrometerreaderFC;
     	// adige.inDams = damsFC;
     	// adige.inTributary = tributaryFC;
     	// adige.inOfftakes = offtakesFC;
     	// adige.inVegetation = vegetationData;
     	adige.pPfafids = "1";
       // adige.pPfafids = "1";
     	adige.fMonpointid = "ID";
     	adige.fPfaff = "pfafstette";
     	adige.fNetelevstart = "elevfirstp";
     	adige.fNetelevend = "elevlastpo";
        adige.pNetNumCali = 9;
        
    //	adige.tTimestep = 60;
		//adige.tStart = startDate;
	//	adige.tEnd = endDate;
		reader_rain.initProcess();
		
      //  adige.process();
		
		adige.pathToDischargeOutput = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/output/discharge_test.csv";
		adige.timeStepMinutes = 60;
		adige.tStartDate = startDate;
		adige.tEndDate = endDate;
		adige.pm = pm;
        adige.process();
		
	     while(reader_rain.doProcess)  {
				reader_rain.nextRecord();
				//System.out.println(reader_rain.tCurrent);
				HashMap<Integer, double[]> rainvalues = reader_rain.outData;
				reader_et.nextRecord();
				HashMap<Integer, double[]> etvalues = reader_et.outData;
	        
				adige.pm = pm;
	            adige.process();
	         
	            
	        
	        }
	   
	   
	        reader_rain.close();
	        reader_et.close();
	       // reader_hydro.close();
	     //   writer.close();
	        System.out.println("ciao");    
	      
	    
	}

	
	
    
    
}



*/