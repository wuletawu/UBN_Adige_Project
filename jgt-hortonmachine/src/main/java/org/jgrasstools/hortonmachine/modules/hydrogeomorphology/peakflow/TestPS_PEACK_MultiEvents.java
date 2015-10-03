package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.peakflow;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.libs.monitor.PrintStreamProgressMonitor;
import org.jgrasstools.hortonmachine.modules.basin.topindex.OmsTopIndex;
import org.jgrasstools.hortonmachine.modules.demmanipulation.pitfiller.OmsPitfiller;
import org.jgrasstools.hortonmachine.modules.geomorphology.draindir.OmsDrainDir;
import org.jgrasstools.hortonmachine.modules.geomorphology.flow.OmsFlowDirections;
import org.jgrasstools.hortonmachine.modules.geomorphology.slope.OmsSlope;
import org.jgrasstools.hortonmachine.modules.network.extractnetwork.OmsExtractNetwork;
import org.jgrasstools.hortonmachine.modules.statistics.kriging.PS;
import org.jgrasstools.hortonmachine.utils.HMTestCase;

/**
 * Test the {@link Tca} module.
 * 
 * @author Giuseppe Formetta ()
 */
public class TestPS_PEACK_MultiEvents extends HMTestCase {

	public void testDream() throws Exception {

		String path_dem = new File(
				"/Users/giuseppeformetta/Desktop/AlbaneGiuseppe/Location/ID4000/cell/elevation_EPSG2154_id4000")
				.getAbsolutePath();

		OmsRasterReader reader_dem = new OmsRasterReader();
		reader_dem.file = path_dem;
		reader_dem.fileNovalue = -9999.0;
		reader_dem.geodataNovalue = Double.NaN;
		reader_dem.process();
		GridCoverage2D dem = reader_dem.outRaster;

		OmsPitfiller pitfiller = new OmsPitfiller();
		pitfiller.inElev = dem;
		pitfiller.pm = pm;
		pitfiller.process();

		GridCoverage2D pitfillerCoverage = pitfiller.outPit;

		OmsFlowDirections flowDirections = new OmsFlowDirections();
		flowDirections.inPit = pitfillerCoverage;
		flowDirections.pm = pm;

		flowDirections.process();

		GridCoverage2D flowCoverage = flowDirections.outFlow;

		OmsDrainDir drainDir = new OmsDrainDir();
		// drainDir.doLad = false;
		drainDir.pLambda = 1;
		drainDir.inPit = pitfillerCoverage;
		drainDir.inFlow = flowCoverage;
		drainDir.pm = pm;

		drainDir.process();

		GridCoverage2D draindirCoverage = drainDir.outFlow;
		GridCoverage2D tcaCoverage = drainDir.outTca;

		OmsSlope slope = new OmsSlope();
		slope.inPit = pitfillerCoverage;
		slope.inFlow = flowCoverage;
		slope.pm = pm;

		slope.process();

		GridCoverage2D slopeCoverage = slope.outSlope;

		OmsTopIndex topindex = new OmsTopIndex();
		topindex.inTca = tcaCoverage;
		topindex.inSlope = slopeCoverage;
		topindex.pm = pm;
		topindex.process();

		GridCoverage2D topindexCoverage = topindex.outTopindex;

		OmsExtractNetwork extractNetwork = new OmsExtractNetwork();
		extractNetwork.pm = pm;
		extractNetwork.inFlow = draindirCoverage;
		extractNetwork.inTca = tcaCoverage;
		extractNetwork.pMode = "tca and slope";
		extractNetwork.pThres = 20000;
		extractNetwork.process();

		GridCoverage2D networkCoverage = extractNetwork.outNet;

		String pathRain = "/Users/giuseppeformetta/Desktop/AlbaneGiuseppe/Events/Rainfall";

		File folder = new File(pathRain);
		File[] listOfFiles = folder.listFiles();
		for (int i = 1; i < listOfFiles.length; i++) {
			System.out.println(listOfFiles[i].getName());
			// vet[i-1]=java.util.Arrays.toString(listOfFiles[i].getName().split(".asc")).toString();
			// System.out.println(vet[i-1]);

		}

		String pathRunoff = "/Users/giuseppeformetta/Desktop/AlbaneGiuseppe/Events/Runoff";

		File folderRunoff = new File(pathRunoff);
		File[] listOfFilesRunoff = folderRunoff.listFiles();
		for (int i = 1; i < listOfFilesRunoff.length; i++) {
			System.out.println(listOfFilesRunoff[i].getName());
			// vet[i-1]=java.util.Arrays.toString(listOfFiles[i].getName().split(".asc")).toString();
			// System.out.println(vet[i-1]);

		}

		// String path = new File(
		//
		// "/Users/giuseppeformetta/Desktop/VAlmazia_GeotopNewAge/locatio_1/val_mazia/cell/basinnumbered1")
		// .getAbsolutePath();
		// RasterReader reader = new RasterReader();
		// reader.file = path;
		// reader.fileNovalue = -9999.0;
		// reader.geodataNovalue = Double.NaN;
		// reader.process();
		// GridCoverage2D rasterbacini = reader.geodata;

		for (int ii = 0; ii < listOfFiles.length; ii++) {

			String pathRainEvent = listOfFiles[ii].toString();
			String pathRunoffEvent = listOfFilesRunoff[ii].toString();
			System.out.println("Event Rainfall= " + pathRainEvent
					+ "  Event Runoff= " + pathRunoffEvent);

			String a = "/Users/giuseppeformetta/Desktop/AlbaneGiuseppe/Events/Par/ParametersEvent_";
			String ppp = a.concat(Integer.toString(ii+1));
			System.out.println(ppp);

			PS d = new PS();
			
			 double[] pmin = { 1.0, 50.0, 0.6, 500.0, 0.1, 0.2 };
			 double[] pmax = { 20, 500.0, 2.0, 10000, 0.9, 1.5 };
			
			 d.ParRange_minn = pmin;
			 d.ParRange_maxn = pmax;
			
			 d.pathToRain = "pathRainEvent";
			 d.pathToRunoff = "pathRunoffEvent";
			 d.ModelName = "PeakFlow";
			 PrintStreamProgressMonitor pm = new PrintStreamProgressMonitor(
			 System.out, System.out);
			 d.pm = pm;
			 d.inNet = networkCoverage;
			 d.inTopIndex = topindexCoverage;
			 d.inDrainDir = draindirCoverage;
			 d.kmax = 15000;
			 d.p = 20;
			 d.parameters = 6;
			
			 d.process();

			double[] vettoprint = d.g_best;

			FileWriter Rstatfile = new FileWriter(ppp);
			PrintWriter errestat = new PrintWriter(Rstatfile);
			for (int j = 0; j < (vettoprint.length); j++) {

				errestat.print(vettoprint[j] + " ");
				System.out.print(vettoprint[j] + " ");

			}
			errestat.println(d.g_best_value+" "+(ii+1));
			System.out.println();

			Rstatfile.close();

		}
	}
}