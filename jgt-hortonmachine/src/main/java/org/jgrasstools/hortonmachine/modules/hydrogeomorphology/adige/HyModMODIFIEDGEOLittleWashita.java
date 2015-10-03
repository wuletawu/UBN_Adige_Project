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
 
package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige;

import java.util.HashMap;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.shapefile.ShapefileFeatureReader;
import org.jgrasstools.gears.io.timedependent.TimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.TimeSeriesIteratorWriter;
import org.jgrasstools.hortonmachine.utils.HMTestCase;

*//**
 * Test ab.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 *//*
public class HyModMODIFIEDGEOLittleWashita extends HMTestCase {
	public static final double doubleNovalue = Double.NaN;

	public void testAb() throws Exception {

		ShapefileFeatureReader hydrometersReader = new ShapefileFeatureReader();
		hydrometersReader.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/hydro4.shp";
		hydrometersReader.readFeatureCollection();
		SimpleFeatureCollection hydrometersFC = hydrometersReader.geodata;

		ShapefileFeatureReader networkReader = new ShapefileFeatureReader();
		networkReader.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/pfaf_tca400_mode2_6.shp";
		networkReader.readFeatureCollection();
		SimpleFeatureCollection networkFC = networkReader.geodata;

		// ShapefileFeatureReader basin5000Reader = new
		// ShapefileFeatureReader();
		// basin5000Reader.file =
		// "/Users/formeppe/Desktop/ORARIO/subbasin50000.shp";
		// basin5000Reader.readFeatureCollection();
		// SimpleFeatureCollection subbasin5000 = basin5000Reader.geodata;

		ShapefileFeatureReader basin5000Reader = new ShapefileFeatureReader();
		basin5000Reader.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/subbasin_tca4000_mode2_6.shp";
		basin5000Reader.readFeatureCollection();
		SimpleFeatureCollection subbasin5000 = basin5000Reader.geodata;

		//
		// /////////////////////////////////////////////////////////////////////////////////////////////////////////
		//
		// /////////////////////////////////////ADIGE////////////////////////////////////////////////////
		//
		// /////////////////////////////////////////////////////////////////////////////////////////////////////////

		String startDate = "2002-01-01 00:00";
		String endDate = "2002-01-02 00:00";
		int timeStepMinutes = 60;

		String fId = "ID";

		TimeSeriesIteratorReader reader_rain = new TimeSeriesIteratorReader();
		reader_rain.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/output/P_DK_interpolated.csv";
		// reader_rain.file =
		// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

		reader_rain.idfield = "ID";
		reader_rain.tStart = startDate;
		reader_rain.tEnd = endDate;
		reader_rain.fileNovalue = "-999.0";
		reader_rain.tTimestep = timeStepMinutes;

		TimeSeriesIteratorWriter writer = new TimeSeriesIteratorWriter();
		writer.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/output/discharge_test.csv";

		writer.tStart = reader_rain.tStart;
		writer.tTimestep = reader_rain.tTimestep;

		TimeSeriesIteratorWriter writer3 = new TimeSeriesIteratorWriter();
		writer3.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/output/discharge_test3.csv";

		writer3.tStart = reader_rain.tStart;
		writer3.tTimestep = reader_rain.tTimestep;

		TimeSeriesIteratorWriter writer2 = new TimeSeriesIteratorWriter();
		writer2.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/output/discharge_test2.csv";

		writer2.tStart = reader_rain.tStart;
		writer2.tTimestep = reader_rain.tTimestep;

		TimeSeriesIteratorReader reader_hydro = new TimeSeriesIteratorReader();
		reader_hydro.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/Hydrometers1995_2013_Novalue.csv";
		reader_hydro.idfield = "ID";
		reader_hydro.tStart = startDate;
		reader_hydro.tEnd = endDate;
		reader_hydro.fileNovalue = "-9999";
		reader_hydro.tTimestep = timeStepMinutes;

		TimeSeriesIteratorReader reader_et = new TimeSeriesIteratorReader();
		reader_et.file = "/Users/administrator/Documents/posina/ISPRAoms3.prj.newage2_etp/data/krigingg_interpolated_etp2.csv";
		reader_et.idfield = "ID";
		reader_et.tStart = startDate;
		reader_et.tEnd = endDate;
		reader_et.fileNovalue = "-999";
		reader_et.tTimestep = timeStepMinutes;

//		HymodInputs hymodInputs = new HymodInputs();
//		hymodInputs.pCmax = 106.55212488441822;
//		hymodInputs.pB = 0.971346077530604;
//		hymodInputs.pAlpha = 0.8084184792349338;
//		hymodInputs.pRs = 89.24261525315367;
//		hymodInputs.pRq = 500.40606747636673;
//		hymodInputs.pQ0 = 0.1;
//
//		hymodInputs.pDorouting = false;
//		hymodInputs.pDoGeom = true;
//		hymodInputs.fIdDistanze = "dtoou";
//		hymodInputs.fIdnetnum = "netnum";
//		hymodInputs.inSubbasinDist = subbasin5000;

		// hymodInputs.betanonlinear=0.9907088317920661;
		// HymodInputs hymodInputs = new HymodInputs();
		// hymodInputs.pCmax =999.9999261305667;
		// hymodInputs.pB =0.5132375402701522 ;
		// hymodInputs.pAlpha =0.23019744642557227;
		// hymodInputs.pRs = 0.015475866372457908;
		// hymodInputs.pRq =0.10436940271552328;
		// hymodInputs.pQ0 = 0.1;
		// hymodInputs.vR=1.49880129653624634;
		// hymodInputs.QR=1.0;
		// hymodInputs.LAMBDA1=0.07472103208709549 ;
		// hymodInputs.LAMBDA2= 0.40687617309471785;
		// hymodInputs.pDorouting=false;
		AdigeModified adige = new AdigeModified();
		adige.pm = pm;
		adige.pCmax = 106.55212488441822;
		adige.pB = 0.971346077530604;
		adige.pAlpha = 0.8084184792349338;
		adige.pRs = 89.24261525315367;
		adige.pRq = 500.40606747636673;
		adige.pQ0 = 0.1;

		adige.pDorouting = false;
		adige.pDoGeom = true;
		adige.pmode=1;
		adige.fIdDistanze = "dtoou";
		adige.fIdnetnum = "netnum";
		adige.inSubbasinDist = subbasin5000;		
		adige.inHillslope = subbasin5000;
		adige.fNetnum = "netnum";
		adige.fBaricenter = "avgZ";
		// adige.fVegetation = "uso_reclas";
		adige.inHydrometers = hydrometersFC;
		// adige.inDams = damsFC;
		// adige.inTributary = tributaryFC;
		// adige.inOfftakes = offtakesFC;
		// adige.inVegetation = vegetationData;
		adige.pPfafids = "1";
		// adige.pPfafids = "1";
		adige.fMonpointid = "ID";
		adige.inNetwork = networkFC;
		adige.fPfaff = "pfafstette";
		adige.fNetelevstart = "elevfirstp";
		adige.fNetelevend = "elevlastpo";

		adige.pRainintensity = -1;
		adige.pRainduration = -1;
		adige.doLog = false;
		adige.tTimestep = timeStepMinutes;
		adige.tStart = startDate;
		adige.tEnd = endDate;
		int indmod = 0;
		adige.pDoReadMeas = true;
		adige.pPathtoMeas = "/Users/administrator/Desktop/test_dischareg";
		adige.pDimMeas = 25;
		adige.pNetNumCali = 1;
		reader_rain.initProcess();
		while (reader_rain.doProcess) {
			reader_rain.nextRecord();
			HashMap<Integer, double[]> rainvalues = reader_rain.outData;
			reader_et.nextRecord();
			HashMap<Integer, double[]> etvalues = reader_et.outData;
			adige.inRain = rainvalues;
			adige.inEtp = etvalues;

			//
			reader_hydro.nextRecord();
			adige.inHydrometerdata = reader_hydro.outData;

			adige.process();

			HashMap<Integer, double[]> outDischarge = adige.outDischarge;
			writer.inData = outDischarge;
			writer.writeNextLine();

			HashMap<Integer, double[]> outvel = adige.outVelocity;
			writer3.inData = outvel;
			writer3.writeNextLine();

			indmod++;
			HashMap<Integer, double[]> outSubDischarge = adige.outSubdischarge;
			writer2.inData = outSubDischarge;
			writer2.writeNextLine();
		}

		reader_rain.close();
		reader_et.close();
		reader_hydro.close();
		writer.close();
		writer2.close();
		writer3.close();
		//
		// /////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ///////////////////////////////////////FINE
		// ADIGE////////////////////////////////////////////////////
		//
		// /////////////////////////////////////////////////////////////////////////////////////////////////////////

		System.out.println("ciao");

	}
}

*/