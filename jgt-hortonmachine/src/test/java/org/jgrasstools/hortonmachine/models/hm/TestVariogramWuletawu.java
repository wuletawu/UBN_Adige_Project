package org.jgrasstools.hortonmachine.models.hm;

import java.util.HashMap;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.jgrasstools.hortonmachine.modules.statistics.kriging.KrigingRagInf012014;
import org.jgrasstools.hortonmachine.modules.statistics.kriging.PS;
import org.jgrasstools.hortonmachine.modules.statistics.kriging.Variogram;
import org.jgrasstools.hortonmachine.utils.HMTestCase;

public class TestVariogramWuletawu extends HMTestCase {
	@SuppressWarnings("nls")
	public void testVariogram() throws Exception {

		String tstart = "2002-01-01 00:00";
		// String tstart = "1980-12-03 00:00";
		String tend = "2004-01-01 00:00";

		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "/Users/administrator/Documents/posina/Temp.prj.newage/data/final_station.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		OmsShapefileFeatureReader pointReader = new OmsShapefileFeatureReader();
		pointReader.file = "/Users/administrator/Documents/posina/Temp.prj.newage/data/centroid_tca400_mode2.shp";
		pointReader.readFeatureCollection();
		SimpleFeatureCollection pointFC = pointReader.geodata;

		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = "/Users/administrator/Documents/posina/Temp.prj.newage/data/T1994_2013.csv";
		reader.idfield = "ID_T";
		reader.tStart = tstart;
		reader.tTimestep = 60;
		reader.tEnd = tend;
		reader.fileNovalue = "-9999";

		reader.initProcess();

		// RasterReader readerr = new RasterReader();
		// readerr.file =
		// "/Users/giuseppeformetta/Desktop/paperKampft/Loc26913/mapset/cell/Aspect300";
		// readerr.fileNovalue = -9999.0;
		// readerr.geodataNovalue = Double.NaN;
		// readerr.process();
		// GridCoverage2D readCoverage = readerr.outRaster;

		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = "/Users/administrator/Documents/posina/Temp.prj.newage/data/INterpolated2002_2004.csv";
		writer.fileNovalue = "-9999.0";
		writer.tStart = tstart;
		writer.tTimestep = 60;
		
		OmsTimeSeriesIteratorWriter writer2 = new OmsTimeSeriesIteratorWriter();
		writer2.file = "/Users/administrator/Documents/posina/Temp.prj.newage/data/INterpolated2002_2004_VAR.csv";
		writer2.fileNovalue = "-9999.0";
		writer2.tStart = tstart;
		writer2.tTimestep = 60;

		String vetmodel[] = new String[] {"gaussian", "linear", "pentaspherical", "bessel", "exponential", "hole" };// ,"pentaspherical"};//;"linear"};//,"bessel","hole", gaussian};

		Variogram Meuse = new Variogram();
		Meuse.pm = pm;
		Meuse.inStations = stationsFC;
		Meuse.fStationsid = "ID_T";
		Meuse.pCutoff = 100000;

		KrigingRagInf012014 kriging = new KrigingRagInf012014();
		kriging.inInterpolate = pointFC;
		kriging.pm = pm;

		kriging.inStations = stationsFC;
		kriging.fStationsid = "ID_T";
		kriging.fInterpolateid = "netnum";
		// kriging.inGridCoverage2D = readCoverage;

		// it doesn't execute the model with log value.
		kriging.doLogarithmic = false;
		/*
		 * Set up the model in order to use the variogram with an explicit
		 * integral scale and variance.
		 */
		// kriging.pVariance = 3.5;
		// kriging.pIntegralscale = new double[]{10000, 10000, 100};

		/*
		 * Set up the model in order to run with a FeatureCollection as point to
		 * interpolated. In this case only 2D.
		 */
		kriging.pMode = 0;
		int iii = 0;
		while (reader.doProcess) {
			reader.nextRecord();
			System.out.println(reader.tCurrent);
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			Meuse.inData = id2ValueMap;
			Meuse.process();
			if (Meuse.outAllEquals == false && Meuse.diversi > 2) {
				// for (int i = 0; i < Meuse.outVar.length; i++) {
				// System.out.println(i + " " + Meuse.outNumPairs[i] + " "
				// + Meuse.outDist[i] + " " + Meuse.outVar[i]);
				// }
				double bestsill = 0;
				double bestnugget = 0;
				double bestraange = 0;
				double minRMSE = 0;
				String bestm = "pippo";
				
				PS a = new PS();

				double maxvar = Meuse.outVar[0];
				double maxdist = Meuse.outDist[0];
				double minvar = Meuse.outVar[0];
				double mindist = Meuse.outDist[0];
				for (int h = 1; h < Meuse.outVar.length; h++) {
					if (Meuse.outVar[h] > maxvar) {
						maxvar = Meuse.outVar[h];
					}
					if (Meuse.outVar[h] < minvar) {
						minvar = Meuse.outVar[h];
					}
					if (Meuse.outDist[h] > maxdist) {
						maxdist = Meuse.outDist[h];
					}
					if (Meuse.outDist[h] < mindist) {
						mindist = Meuse.outDist[h];
					}
				}

				for (int ooo = 0; ooo < vetmodel.length; ooo++) {

					a.inVar = Meuse.outVar;
					a.inDistance = Meuse.outDist;
					a.inmodelname = vetmodel[ooo];
					a.ModelName = "Vgm";
					a.parameters = 3;
					a.kmax = 10000;
					a.threshold = 1e-8;
					a.p = 15;
					a.inNp = Meuse.outNumPairs;
					a.ParRange_minn = new double[] { 0.0, 0.0, mindist };
					a.ParRange_maxn = new double[] { 0.6, 6.8, 110000 };
					a.process();
					// System.out.println(a.g_best_value + " nugget="
					// + a.g_best[0] + " sill=" + a.g_best[1] + " range="
					// + a.g_best[2]);
					if (ooo == 0) {
						bestm = vetmodel[ooo];
						bestsill = a.g_best[1];
						bestraange = a.g_best[2];
						bestnugget = a.g_best[0];
						minRMSE = a.g_best_value;

					} else {
						if (a.g_best_value < minRMSE) {
							bestm = vetmodel[ooo];
							bestsill = a.g_best[1];
							bestraange = a.g_best[2];
							bestnugget = a.g_best[0];
							minRMSE = a.g_best_value;
						}

					}

				}
				kriging.doDetrended = true;
				kriging.plocalTrend = false;
				kriging.fPointZ = "avgZ";
				 kriging.thresholdCorrelation = -0.9;
				// kriging.maxdist = 100;
				 kriging.inNumCloserStations=5;
				// kriging.fStationsZ = "int_1";
				kriging.defaultVariogramMode = 1;
				kriging.pA = bestraange;
				kriging.pNug = bestnugget;
				kriging.pS = bestsill;

				// kriging.doIncludezero = true;
				// kriging.constrainGT0 = true;
				kriging.inData = id2ValueMap;
				kriging.pSemivariogramType = bestm;
				kriging.executeKriging();
				System.out
						.println("HO EFFETTUATO IL KRIGING CON QUESTI PARAMETRI:"
								+ " nugget="
								+ kriging.pNug
								+ " range="
								+ kriging.pA
								+ " sill="
								+ kriging.pS
								+ " model=" + kriging.pSemivariogramType);
				// RasterWriter writer = new RasterWriter();
				// writer.inRaster = kriging.outGrid;
				// writer.file =
				// "/Users/giuseppeformetta/Desktop/paperKampft/Loc26913/mapset/cell/T_"
				// + iii;
				// writer.process();
				// iii++;

				writer.inData = kriging.outData;
				writer.writeNextLine();
				
				writer2.inData = kriging.outVarData;
				writer2.writeNextLine();

			} else {
				kriging.doDetrended = true;
				kriging.thresholdCorrelation = -0.9;
				// kriging.plocalTrend=true;

				 kriging.inNumCloserStations = 5;
				// kriging.maxdist = 1000;
				// kriging.doIncludezero=false;
				// kriging.constrainGT0 = true;
				kriging.fStationsZ = "elv";
				kriging.defaultVariogramMode = 1;
				kriging.pA = 100000.0;
				kriging.pNug = 0.1;
				kriging.pS = 0.1;
				kriging.inData = id2ValueMap;
				kriging.pSemivariogramType = "linear";
				kriging.executeKriging();
				System.out
						.println("HO EFFETTUATO IL KRIGING CON QUESTI PARAMETRI:"
								+ " nugget="
								+ kriging.pNug
								+ " range="
								+ kriging.pA
								+ " sill="
								+ kriging.pS
								+ " model=" + kriging.pSemivariogramType);
				// RasterWriter writer = new RasterWriter();
				// writer.inRaster = kriging.outGrid;
				// writer.file =
				// "/Users/giuseppeformetta/Desktop/paperKampft/Loc26913/mapset/cell/T_"
				// + iii;
				// writer.process();

				writer.inData = kriging.outData;
				writer.writeNextLine();
				iii++;
				
				writer2.inData = kriging.outVarData;
				writer2.writeNextLine();
				iii++;
			}

		}
		System.out.println("ciao");
		writer.close();
		writer2.close();

	}
}
