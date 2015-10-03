package org.jgrasstools.hortonmachine.modules.statistics.kriging;

import static org.jgrasstools.gears.libs.modules.JGTConstants.doubleNovalue;

import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Set;

import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;
import oms3.annotations.Unit;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timeseries.OmsTimeSeriesReader;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.monitor.DummyProgressMonitor;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.hortonmachine.modules.basin.rescaleddistance.OmsRescaledDistance;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.AdigeModified;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.peakflow.OmsPeakflow;
import org.joda.time.DateTime;
//import corejava.Format;
//import sun.tools.asm.Cover;
//import umontreal.iro.lecuyer.probdist.*;
//import umontreal.iro.lecuyer.stat.*;

public class PS extends JGTModel {

	@Description("The progress monitor.")
	@In
	public String START_DATE = null;

	@Description("The map of the skyview factor.")
	@In
	public GridCoverage2D inskyview = null;

	@Description("The vector of the measurement point, containing the position of the stations.")
	@In
	public SimpleFeatureCollection inStations = null;

	@Description("The field of the vector of stations, defining the id.")
	@In
	public String fStationsid = null;

	@Description("The field of the vector of stations, defining the id.")
	@In
	public GridCoverage2D inNet = null;

	@Description("The field of the vector of stations, defining the id.")
	@In
	public GridCoverage2D inTopIndex = null;

	@Description("The field of the vector of stations, defining the id.")
	@In
	public GridCoverage2D inDrainDir = null;

	@Description("The progress monitor.")
	@In
	public String END_DATE = null;

	@Description("The map of the dem.")
	@In
	public GridCoverage2D dem = null;

	@Description("The map of the insolation for February.")
	@In
	public GridCoverage2D insolationFebruary = null;

	@Description("The map of the insolation for March.")
	@In
	public GridCoverage2D insolationMarch = null;

	@Description("The map of the insolation for March.")
	@In
	public double[] inVar = null;

	@Description("The map of the insolation for March.")
	@In
	public double[] inNp = null;

	@Description("The map of the insolation for March.")
	@In
	public double[] inDistance = null;
	@Description("The map of the insolation for March.")
	@In
	public String inmodelname = null;

	@Description("The map of the insolation for April.")
	@In
	public GridCoverage2D insolationApril = null;

	@Description("The map of the insolation for May.")
	@In
	public GridCoverage2D insolationMay = null;

	@Description("The map of the insolation for June.")
	@In
	public GridCoverage2D insolationJune = null;
	// ////////////////////////////////////////////
	@Description("The map of the skyview factor.")
	@In
	public GridCoverage2D inSkyview = null;
	@Description("The first day of the simulation.")
	@In
	public String tStartDate = null;
	@Description("The time step in minutes of the measurement.")
	@In
	public int inTimestep;
	@Description("The last day of the simulation.")
	@In
	public String tEndDate = null;

	@Description("The map of the dem.")
	@In
	public GridCoverage2D inDem = null;

	@Description("The first day of the simulation.")
	@In
	public String time = null;

	@Description("The Tmelt.")
	@In
	@Unit("C")
	public double pTmelt;
	@Description("Path to the temperature file in input.")
	@In
	public String pPathToRain = null;

	@Description("Path to the humidity file in input.")
	@In
	public String pPathToTemp = null;
	@Description("Path to the temperature file in input.")
	@In
	public String pPathToCondIniI = null;

	@Description("Path to the humidity file in input.")
	@In
	public String pPathToCondIniL = null;

	@Description("The pr of lmax.")
	@In
	@Unit("C")
	public double pR;
	@Description("The path to direct output.")
	@In
	public String pPathSWE;
	@Description("The path to direct output.")
	@In
	public String pPathMelting;
	@Description("The map of the insolation for February.")
	@In
	public GridCoverage2D inInsFeb = null;

	@Description("The map of the insolation for March.")
	@In
	public GridCoverage2D inInsMar = null;

	@Description("The map of the insolation for April.")
	@In
	public GridCoverage2D inInsApr = null;

	@Description("The map of the insolation for May.")
	@In
	public GridCoverage2D inInsMay = null;

	@Description("The map of the insolation for June.")
	@In
	public GridCoverage2D inInsJun = null;
	@Description("The map of the insolation for June.")
	@In
	public double threshold = 0.001;

	@Description("Do you want raster map as output.")
	@In
	public boolean doRaster = true;

	@Description("The temperature hashmap.")
	@In
	@Unit(" C ")
	public HashMap<Integer, double[]> inTemp;

	@Description("The rainfall.")
	@In
	@Unit("mm")
	public HashMap<Integer, double[]> inRainfall;

	@Description("The inital condition I.")
	@In
	@Unit("mm")
	public HashMap<Integer, double[]> inCondIniI;
	@Description("The inital condition L")
	@In
	@Unit("mm")
	public HashMap<Integer, double[]> inCondIniL;
	// ////////////////////////////////////////////////////

	@Description("The progress monitor.")
	@In
	public IJGTProgressMonitor pm = new DummyProgressMonitor();

	@Description("")
	@In
	public int Warmup;

	@Description("Potential ETP vector.")
	@In
	public double[] Extra_PET;

	@Description("Precipitations vector.")
	@In
	public double[] Extra_Precip;

	@Description("The limit of the time series considered in computation.")
	@In
	public double Extra_MaxT = doubleNovalue;

	@Description("The network data.")
	@In
	public SimpleFeatureCollection net;

	@Description("The network data.")
	@In
	public SimpleFeatureCollection hydrome;

	@Description("The network data.")
	@In
	public SimpleFeatureCollection basin;

	@Description("Define the measured streamflow data.")
	@In
	public double[] Measurement_MeasData;

	@Description("Area of the basin.")
	@In
	public double Area;

	@Description("The first value of mesured discharge")
	@In
	@Unit("m3/s")
	public double Q0 = doubleNovalue;

	@Description("kmax")
	@In
	public int kmax;

	@Description("Rain input path")
	@In
	public String pathToRain;

	@Description("Runoff input path")
	@In
	public String pathToRunoff;

	@Description("ModelName")
	@In
	public String ModelName = null;

	@Description("Give the parameter ranges (minimum values).")
	@In
	public double[] ParRange_minn;

	@Description("Give the parameter ranges (maximum values).")
	@In
	public double[] ParRange_maxn;

	@Description("p")
	@In
	public int p;
	@Description("parameters")
	@In
	public int parameters;

	@Description("gbest")
	@Out
	public double outOpt;

	@Description("optimal parameters")
	@Out
	public double[] outOptSet;
	@Description("optimal parameters")
	@Out
	public double[] observedvect;

	public int int_Extra_MaxT;
	public double[] costvect;
	public double[][] p_best;
	public double[] g_best;
	public double g_best_value;
	public double[] p_best_value;
	public double[] hist_g_best_value;
	public double observed[];
	public double modelled[];

	@Execute
	public void process() throws Exception {
		int_Extra_MaxT = (int) Extra_MaxT;

		double xX[][] = new double[parameters][p];

		xX = UNIFORM(ParRange_minn, ParRange_maxn, p);
		double x[][] = ReflectBounds(xX, ParRange_maxn, ParRange_minn);
		// for (int i = 0; i < x.length; i++) {
		// for (int j = 0; j < x[0].length;j++){
		// System.out.print(x[i][j]+ " ");
		// }
		// System.out.println();
		//
		// }
		double VelRange_minn[] = new double[ParRange_minn.length];

		double VelRange_max[] = new double[ParRange_minn.length];
		// ipotizzo che la velocita iniziale sia inizializzata compresa tra
		// 1/100 dei valori max e min dei parametri
		for (int i = 0; i < ParRange_maxn.length; i++) {
			VelRange_minn[i] = ParRange_minn[i];
			VelRange_max[i] = ParRange_maxn[i];
		}

		double vel[][] = new double[parameters][p];
		vel = UNIFORM(VelRange_minn, VelRange_max, p);

		// calculate the cost of each particle
		double[] costvectold = ComputeCostFunction(x);

		int kkk = 0;
		p_best = x;
		p_best_value = costvectold;
		double min = Math.abs(costvectold[0]);
		int posmin = 0;
		for (int i = 1; i < costvectold.length; i++) {
			if (Math.abs(costvectold[i]) < min) {
				min = Math.abs(costvectold[i]);
				posmin = i;
			}
		}
		g_best = new double[x[0].length];
		g_best_value = min;
		for (int i = 0; i < x[0].length; i++) {
			g_best[i] = x[posmin][i];
		}

		hist_g_best_value = new double[kmax];
		boolean fermati = false;
		while (kkk < (kmax - 1) || !fermati) {
			double[][] x_old = x;
			double[][] velnew = Compute_velocity(x_old, vel);
			vel = velnew;
			x = Compute_particle(x_old, velnew);
			costvect = ComputeCostFunction(x);
			p_best = Compute_pBest(x, costvect);
			g_best = Compute_gBest(x, costvect);

			hist_g_best_value[kkk] = g_best_value;
			for (int jjj = 0; jjj < g_best.length; jjj++) {

				System.out.print(g_best[jjj] + "  ");

			}
			System.out.println("g_best_value= " + g_best_value);

			if (kkk > 1000) {
				int sum = 0;
				for (int c = 0; c < 550; c++) {
					if (Math.abs(hist_g_best_value[kkk - c]
							- hist_g_best_value[kkk - c - 1]) < threshold) {
						sum = sum + 1;
					}
				}

				if (sum > 30) {
					fermati = true;
					break;
				}
			}
			if (kkk > kmax - 2) {
				break;
			}
			// aggiunto per NASH ALBAN
			// if (ModelName.equals("PeakFlow")) {
			// if (g_best_value < 0.05) {
			// break;
			// }
			// }
			costvectold = costvect;
			// System.out.println("ite="+kkk);
			// System.out.println("fermati="+fermati);
			kkk++;

		}
		// for (int i = 0; i < g_best.length; i++) {
		// System.out.println(g_best[i]);
		//
		// }
		// System.out.println(g_best_value);

		outOpt = g_best_value;
		outOptSet = g_best;
	}

	public boolean StoppingCondition(double vett[], double s) {

		boolean result = false;
		int cnt = 0;
		for (int i = 0; i < vett.length; i++) {
			if (vett[i] < s) {
				cnt += 1;
			}
		}
		if (cnt == vett.length) {
			result = true;
		}

		return result;
	}

	public double[][] UNIFORM(double[] xmin, double[] xmax, int nsample) {
		// Latin Hypercube sampling
		// double[][] LHSresult=new double [1][1];
		int nvar = xmin.length;
		double[][] s = new double[nsample][nvar];

		double[][] ran = new double[nsample][nvar];
		for (int row = 0; row < ran.length; row++) {

			for (int col = 0; col < ran[0].length; col++) {

				s[row][col] = (xmax[col] - xmin[col]) * Math.random()
						+ xmin[col];
			}
		}

		return s;

	}

	public double[] ComputeCostFunction(double xx[][]) throws Exception {
		double[] res = new double[xx.length];

		//

		// ///////////////////////////////////////////////////////////////////////////////////////////////////////
		// ///////////////////////////////////////F I N E A D I G
		// E//////////////////////////////////////////////////
		// ///////////////////////////////////////////////////////////////////////////////////////////////////////
		/*if (ModelName.equals("AdigeLW")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				String startDate = "2002-01-01 00:00";
				String endDate = "2004-12-31 00:00";
				int timeStepMinutes = 60 * 24;

				OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
				reader_rain.file = "/Users/giuseppeformetta/Desktop/LittleWshita/hourlydataset/Rain2002_2008_Daily.csv";
				// reader_rain.file =
				// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

				reader_rain.idfield = "ID";
				reader_rain.tStart = startDate;
				reader_rain.tEnd = endDate;
				reader_rain.fileNovalue = "-9999.0";
				reader_rain.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
				reader_hydro.file = "/Users/giuseppeformetta/Desktop/LittleWshita/hourlydataset/Hydrometers_DAily2002_2008.csv";
				reader_hydro.idfield = "RETE_ID";
				reader_hydro.tStart = startDate;
				reader_hydro.tEnd = endDate;
				reader_hydro.fileNovalue = "-9999";
				reader_hydro.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
				reader_et.file = "/Users/giuseppeformetta/Desktop/LittleWshita/hourlydataset/Etp2002_2008_Daily.csv";
				reader_et.idfield = "ID";
				reader_et.tStart = startDate;
				reader_et.tEnd = endDate;
				reader_et.fileNovalue = "-9999";
				reader_et.tTimestep = timeStepMinutes;

				int ii = numpart;

				AdigeModified adige = new AdigeModified();
				adige.pCmax = xx[ii][0];
				adige.pB = xx[ii][1];
				adige.pAlpha = xx[ii][2];
				adige.pRs = xx[ii][3];
				adige.pRq = xx[ii][4];
				adige.pMetp = xx[ii][5];
				adige.pMrain = xx[ii][6];
				// adige.pQr = 1.0;
				// adige.vR = xx[ii][7];
				// adige.pLambda1 = xx[ii][8];
				// adige.pLambda2 = xx[ii][9];
				//adige.pNonlinExp = 1.0;
				adige.pQ0 = 0.0;
				//adige.doNonlinearRes = false;
				adige.pDoReadMeas = true;
				adige.pPathtoMeas = "/Users/giuseppeformetta/Desktop/LittleWshita/hourlydataset/LWDAILYOK";
				adige.pDimMeas = 1095;
				adige.pBetanolinear = 1.0;
				adige.pmode = 1;
				adige.pDorouting = false;
				adige.pDoGeom = false;
				adige.fIdDistanze = "dtoou";
				adige.fIdnetnum = "netnum";
				adige.inSubbasinDist = basin;
				adige.pm = pm;
				adige.inHillslope = basin;
				adige.fNetnum = "netnum";
				adige.fBaricenter = "avgZ";
				// adige.fVegetation = "uso_reclas";
				adige.inHydrometers = hydrome;
				// adige.inDams = damsFC;
				// adige.inTributary = tributaryFC;
				// adige.inOfftakes = offtakesFC;
				// adige.inVegetation = vegetationData;
				adige.pPfafids = "3";
				// adige.pPfafids = "1";
				adige.fMonpointid = "RETE_ID";
				adige.inNetwork = net;
				adige.fPfaff = "pfafstette";
				adige.fNetelevstart = "elevfirstp";
				adige.fNetelevend = "elevlastpo";
				// adige.tTimestep = 24;
				adige.pRainintensity = -1;
				adige.pRainduration = -1;
				adige.doLog = false;

				adige.tTimestep = timeStepMinutes;
				adige.tStart = startDate;
				adige.tEnd = endDate;
				adige.pNetNumCali = 2;
				reader_rain.initProcess();
				while (reader_rain.doProcess) {

					reader_rain.nextRecord();
					HashMap<Integer, double[]> rainvalues = reader_rain.outData;
					reader_et.nextRecord();
					HashMap<Integer, double[]> etvalues = reader_et.outData;
					adige.inRain = rainvalues;
					adige.inEtp = etvalues;

					reader_hydro.nextRecord();
					adige.inHydrometerdata = reader_hydro.outData;

					adige.process();

				}

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();
				double[] modellati = new double[730];
				double[] osservati = new double[730];
				for (int i = 365; i < 1095; i++) {
					modellati[i - 365] = adige.outSimulated[i];
					osservati[i - 365] = adige.outMeasured[i];
				}
				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = osservati;
				o.simulated = modellati;
				o.process();
				// va benissimo
				res[numpart] = 1 - (o.outKGE);
				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				System.out.println(numpart + "    " + res[numpart]);

			}

		}*/
		
		if (ModelName.equals("AdigeLW20062008")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				String startDate = "2007-01-01 00:00";
				String endDate = "2008-12-31 00:00";
				int timeStepMinutes = 60;

				OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
				reader_rain.file = "/Users/giuseppeformetta/Desktop/LittleWshita/Rain2002_2008.csv";
				// reader_rain.file =
				// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

				reader_rain.idfield = "ID";
				reader_rain.tStart = startDate;
				reader_rain.tEnd = endDate;
				reader_rain.fileNovalue = "-9999.0";
				reader_rain.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
				reader_hydro.file = "/Users/giuseppeformetta/Desktop/IS/ISPRAoms3.prj.newage/data/Hydrometers2002_2008.csv";
				reader_hydro.idfield = "RETE_ID";
				reader_hydro.tStart = startDate;
				reader_hydro.tEnd = endDate;
				reader_hydro.fileNovalue = "-9999";
				reader_hydro.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
				reader_et.file = "/Users/giuseppeformetta/Desktop/IS/ISPRAoms3.prj.newage/data/Etp2002_2008.csv";
				reader_et.idfield = "ID";
				reader_et.tStart = startDate;
				reader_et.tEnd = endDate;
				reader_et.fileNovalue = "-9999";
				reader_et.tTimestep = timeStepMinutes;

				int ii = numpart;

				AdigeModified adige = new AdigeModified();
				adige.pCmax = xx[ii][0];
				adige.pB = xx[ii][1];
				adige.pAlpha = xx[ii][2];
				adige.pRs = xx[ii][3];
				adige.pRq = xx[ii][4];
				adige.pMetp = xx[ii][5];
				adige.pMrain = xx[ii][6];

				adige.pQ0 = 0.02;
			//	adige.doNonlinearRes = false;
				adige.pDoReadMeas = true;
				adige.pPathtoMeas = "/Users/giuseppeformetta/Desktop/LittleWshita/portate_staz_07327550JAVA20062008";
				adige.pDimMeas = 1*24*365;
				adige.pBetanolinear = 1.0;
				adige.pmode = 1;
				adige.pDorouting = false;
				adige.pDoGeom = false;
				adige.fIdDistanze = "dtoou";
				adige.fIdnetnum = "netnum";
				adige.inSubbasinDist = basin;
				adige.pm = pm;
				adige.inHillslope = basin;
				adige.fNetnum = "netnum";
				adige.fBaricenter = "avgZ";
				// adige.fVegetation = "uso_reclas";
				adige.inHydrometers = hydrome;
				// adige.inDams = damsFC;
				// adige.inTributary = tributaryFC;
				// adige.inOfftakes = offtakesFC;
				// adige.inVegetation = vegetationData;
				adige.pPfafids = "3";
				// adige.pPfafids = "1";
				adige.fMonpointid = "RETE_ID";
				adige.inNetwork = net;
				adige.fPfaff = "pfafstette";
				adige.fNetelevstart = "elevfirstp";
				adige.fNetelevend = "elevlastpo";

				adige.pRainintensity = -1;
				adige.pRainduration = -1;
				adige.doLog = false;

				adige.tTimestep = timeStepMinutes;
				adige.tStart = startDate;
				adige.tEnd = endDate;
				adige.pNetNumCali = 2;
				reader_rain.initProcess();
				while (reader_rain.doProcess) {

					reader_rain.nextRecord();
					HashMap<Integer, double[]> rainvalues = reader_rain.outData;
					reader_et.nextRecord();
					HashMap<Integer, double[]> etvalues = reader_et.outData;
					adige.inRain = rainvalues;
					adige.inEtp = etvalues;

					reader_hydro.nextRecord();
					adige.inHydrometerdata = reader_hydro.outData;

					adige.process();

				}

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = adige.outMeasured;
				o.simulated = adige.outSimulated;
				o.process();
	
				// va benissimo
				res[numpart] = 1 - (o.outKGE);

				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				System.out.println(numpart + "    " + res[numpart]);

			}

		}
		if (ModelName.equals("AdigeLW20062008R")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				String startDate = "2007-01-01 00:00";
				String endDate = "2008-12-31 00:00";
				int timeStepMinutes = 60;

				OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
				reader_rain.file = "/Users/giuseppeformetta/Desktop/LittleWshita/Rain2002_2008.csv";
				// reader_rain.file =
				// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

				reader_rain.idfield = "ID";
				reader_rain.tStart = startDate;
				reader_rain.tEnd = endDate;
				reader_rain.fileNovalue = "-9999.0";
				reader_rain.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
				reader_hydro.file = "/Users/giuseppeformetta/Desktop/IS/ISPRAoms3.prj.newage/data/Hydrometers2002_2008.csv";
				reader_hydro.idfield = "RETE_ID";
				reader_hydro.tStart = startDate;
				reader_hydro.tEnd = endDate;
				reader_hydro.fileNovalue = "-9999";
				reader_hydro.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
				reader_et.file = "/Users/giuseppeformetta/Desktop/IS/ISPRAoms3.prj.newage/data/Etp2002_2008.csv";
				reader_et.idfield = "ID";
				reader_et.tStart = startDate;
				reader_et.tEnd = endDate;
				reader_et.fileNovalue = "-9999";
				reader_et.tTimestep = timeStepMinutes;

				int ii = numpart;

				AdigeModified adige = new AdigeModified();
				adige.pCmax = xx[ii][0];
				adige.pB = xx[ii][1];
				adige.pAlpha = xx[ii][2];
				adige.pRs = xx[ii][3];
				adige.pRq = xx[ii][4];
				adige.pMetp = xx[ii][5];
				adige.pMrain = xx[ii][6];
				adige.vR = xx[ii][7];
				;
				adige.pQr = 1.0;
				adige.pLambda1 = xx[ii][8];
				;
				adige.pLambda2 = xx[ii][9];
				;
				adige.pQ0 = 0.1;
				//adige.doNonlinearRes = false;
				adige.pDoReadMeas = true;
				adige.pPathtoMeas = "/Users/giuseppeformetta/Desktop/LittleWshita/portate_staz_07327550JAVA20062008";
				adige.pDimMeas = 2*24*365;
				adige.pBetanolinear = 1.0;
				adige.pmode = 1;
				adige.pDorouting = true;
				adige.pDoGeom = false;
				adige.fIdDistanze = "dtoou";
				adige.fIdnetnum = "netnum";
				adige.inSubbasinDist = basin;
				adige.pm = pm;
				adige.inHillslope = basin;
				adige.fNetnum = "netnum";
				adige.fBaricenter = "avgZ";
				// adige.fVegetation = "uso_reclas";
				adige.inHydrometers = hydrome;
				// adige.inDams = damsFC;
				// adige.inTributary = tributaryFC;
				// adige.inOfftakes = offtakesFC;
				// adige.inVegetation = vegetationData;
				adige.pPfafids = "3";
				// adige.pPfafids = "1";
				adige.fMonpointid = "RETE_ID";
				adige.inNetwork = net;
				adige.fPfaff = "pfafstette";
				adige.fNetelevstart = "elevfirstp";
				adige.fNetelevend = "elevlastpo";

				adige.pRainintensity = -1;
				adige.pRainduration = -1;
				adige.doLog = false;

				adige.tTimestep = timeStepMinutes;
				adige.tStart = startDate;
				adige.tEnd = endDate;
				adige.pNetNumCali = 2;
				reader_rain.initProcess();
				while (reader_rain.doProcess) {

					reader_rain.nextRecord();
					HashMap<Integer, double[]> rainvalues = reader_rain.outData;
					reader_et.nextRecord();
					HashMap<Integer, double[]> etvalues = reader_et.outData;
					adige.inRain = rainvalues;
					adige.inEtp = etvalues;

					reader_hydro.nextRecord();
					adige.inHydrometerdata = reader_hydro.outData;

					adige.process();

				}

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = adige.outMeasured;
				o.simulated = adige.outSimulated;
				o.process();
	
				// va benissimo
				res[numpart] = 1 - (o.outKGE);

				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				System.out.println(numpart + "    " + res[numpart]);

			}

		}
		
		
		
		if (ModelName.equals("AdigeCobb")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				String startDate = "2006-01-01 00:00";
				String endDate = "2007-01-01 00:00";
				int timeStepMinutes = 60;

				OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
				// reader_rain.file =
				// "/Users/giuseppeformetta/Desktop/COBB/totale_pioggeidw_modif2.csv";
				reader_rain.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/KRIGING_RAIN_DISTRIBUTED.csv";
				// reader_rain.file =
				// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

				reader_rain.idfield = "ID";
				reader_rain.tStart = startDate;
				reader_rain.tEnd = endDate;
				reader_rain.fileNovalue = "-9999.0";
				reader_rain.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
				reader_hydro.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/idrometro.csv";
				reader_hydro.idfield = "Rete_id";
				reader_hydro.tStart = startDate;
				reader_hydro.tEnd = endDate;
				reader_hydro.fileNovalue = "-9999";
				reader_hydro.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
				reader_et.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/totale_ETP_modif.csv";
				reader_et.idfield = "ID";
				reader_et.tStart = startDate;
				reader_et.tEnd = endDate;
				reader_et.fileNovalue = "-9999";
				reader_et.tTimestep = timeStepMinutes;

				int ii = numpart;

				AdigeModified adige = new AdigeModified();
				adige.pCmax = xx[ii][0];
				adige.pB = xx[ii][1];
				adige.pAlpha = xx[ii][2];
				adige.pRs = xx[ii][3];
				adige.pRq = xx[ii][4];
				adige.pMetp = xx[ii][5];
				adige.pMrain = xx[ii][6];
				// adige.pNonlinExp = xx[ii][7];
				adige.pQ0 = 0.2;
				//adige.doNonlinearRes = false;
				adige.pDoReadMeas = true;
				adige.pPathtoMeas = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/MeasuredDischarhe_01012006.csv";
				adige.pDimMeas = 8760;
				adige.pBetanolinear = 1.0;
				adige.pmode = 1;
				adige.pDorouting = false;
				adige.pDoGeom = false;
				adige.fIdDistanze = "dToO";
				adige.fIdnetnum = "netnum";
				adige.inSubbasinDist = basin;
				adige.pm = pm;
				adige.inHillslope = basin;
				adige.fNetnum = "netnum";
				adige.fBaricenter = "avgZ";
				// adige.fVegetation = "uso_reclas";
				adige.inHydrometers = hydrome;
				// adige.inDams = damsFC;
				// adige.inTributary = tributaryFC;
				// adige.inOfftakes = offtakesFC;
				// adige.inVegetation = vegetationData;
				adige.pPfafids = "3";
				// adige.pPfafids = "1";
				adige.fMonpointid = "RETE_ID";
				adige.inNetwork = net;
				adige.fPfaff = "pfafstette";
				adige.fNetelevstart = "elevfirstp";
				adige.fNetelevend = "elevlastpo";

				adige.pRainintensity = -1;
				adige.pRainduration = -1;
				adige.doLog = false;

				adige.tTimestep = timeStepMinutes;
				adige.tStart = startDate;
				adige.tEnd = endDate;
				adige.pNetNumCali = 4;
				reader_rain.initProcess();
				while (reader_rain.doProcess) {

					reader_rain.nextRecord();
					HashMap<Integer, double[]> rainvalues = reader_rain.outData;
					reader_et.nextRecord();
					HashMap<Integer, double[]> etvalues = reader_et.outData;
					adige.inRain = rainvalues;
					adige.inEtp = etvalues;

					reader_hydro.nextRecord();
					adige.inHydrometerdata = reader_hydro.outData;

					adige.process();

				}

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = adige.outMeasured;
				o.simulated = adige.outSimulated;
				o.process();

				// va benissimo
				res[numpart] = 1 - (o.outKGE);
				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				// for(int iii=0;iii<adige.outMeasured.length;iii++){
				// System.out.println( iii+"  "+adige.outMeasured[iii]+ "  "+
				// adige.outSimulated[iii]);
				// }
				System.out.println(numpart + xx[ii][0] + " " + xx[ii][1] + " "
						+ xx[ii][2] + " " + xx[ii][3] + " " + xx[ii][4] + " "
						+ xx[ii][5] + " " + xx[ii][6] + " " + res[numpart]);

			}

		}
		if (ModelName.equals("AdigeLWPDM")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				String startDate = "2002-01-01 00:00";
				String endDate = "2004-12-31 00:00";
				int timeStepMinutes = 60 * 24;
				int ii = numpart;
				OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
				reader_rain.file = "/Users/giuseppeformetta/Desktop/LittleWshita/hourlydataset/Rain2002_2008_Daily.csv";
				// reader_rain.file =
				// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

				reader_rain.idfield = "ID";
				reader_rain.tStart = startDate;
				reader_rain.tEnd = endDate;
				reader_rain.fileNovalue = "-9999.0";
				reader_rain.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
				reader_hydro.file = "/Users/giuseppeformetta/Desktop/LittleWshita/hourlydataset/Hydrometers_DAily2002_2008.csv";
				reader_hydro.idfield = "RETE_ID";
				reader_hydro.tStart = startDate;
				reader_hydro.tEnd = endDate;
				reader_hydro.fileNovalue = "-9999";
				reader_hydro.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
				reader_et.file = "/Users/giuseppeformetta/Desktop/LittleWshita/hourlydataset/Etp2002_2008_Daily.csv";
				reader_et.idfield = "ID";
				reader_et.tStart = startDate;
				reader_et.tEnd = endDate;
				reader_et.fileNovalue = "-9999";
				reader_et.tTimestep = timeStepMinutes;

				AdigeModified adige = new AdigeModified();
				adige.pCmax = xx[ii][0];
				adige.pB = xx[ii][1];
				adige.pAlpha = xx[ii][2];
				adige.pRs = xx[ii][3];
				adige.pRq = xx[ii][4];
				adige.pMetp = xx[ii][5];
				adige.pMrain = xx[ii][6];
				// adige.pQr = 1.0;
				// adige.vR = xx[ii][7];
				// adige.pLambda1 = xx[ii][8];
				// adige.pLambda2 = xx[ii][9];
			//	adige.pNonlinExp = 1.0;
				adige.pQ0 = 0.0;
			//	adige.doNonlinearRes = false;
				adige.pDoReadMeas = true;
				adige.pPathtoMeas = "/Users/giuseppeformetta/Desktop/LittleWshita/hourlydataset/LWDAILYOK";
				adige.pDimMeas = 1095;
				adige.pBetanolinear = 1.0;
				adige.pmode = 2;
				adige.pDorouting = false;
				adige.pDoGeom = false;
				adige.fIdDistanze = "dtoou";
				adige.fIdnetnum = "netnum";
				adige.inSubbasinDist = basin;
				adige.pm = pm;
				adige.inHillslope = basin;
				adige.fNetnum = "netnum";
				adige.fBaricenter = "avgZ";
				// adige.fVegetation = "uso_reclas";
				adige.inHydrometers = hydrome;
				// adige.inDams = damsFC;
				// adige.inTributary = tributaryFC;
				// adige.inOfftakes = offtakesFC;
				// adige.inVegetation = vegetationData;
				adige.pPfafids = "3";
				// adige.pPfafids = "1";
				adige.fMonpointid = "RETE_ID";
				adige.inNetwork = net;
				adige.fPfaff = "pfafstette";
				adige.fNetelevstart = "elevfirstp";
				adige.fNetelevend = "elevlastpo";
				// adige.tTimestep = 24;
				adige.pRainintensity = -1;
				adige.pRainduration = -1;
				adige.doLog = false;

				adige.tTimestep = timeStepMinutes;
				adige.tStart = startDate;
				adige.tEnd = endDate;
				adige.pNetNumCali = 2;
				reader_rain.initProcess();
				while (reader_rain.doProcess) {

					reader_rain.nextRecord();
					HashMap<Integer, double[]> rainvalues = reader_rain.outData;
					reader_et.nextRecord();
					HashMap<Integer, double[]> etvalues = reader_et.outData;
					adige.inRain = rainvalues;
					adige.inEtp = etvalues;

					reader_hydro.nextRecord();
					adige.inHydrometerdata = reader_hydro.outData;

					adige.process();

				}

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();
				// double[] modellati = new double[730];
				// double[] osservati = new double[730];
				// for (int i = 365; i < 1095; i++) {
				// modellati[i - 365] = adige.outSimulated[i];
				// osservati[i - 365] = adige.outMeasured[i];
				// }
				//

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = adige.outMeasured;
				o.simulated = adige.outSimulated;
				o.process();
				// va benissimo
				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				res[numpart] = 1 - (o.outObFunNash);
				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				// for(int iii=0;iii<adige.outMeasured.length;iii++){
				// System.out.println( iii+"  "+adige.outMeasured[iii]+ "  "+
				// adige.outSimulated[iii]);
				// }
				System.out.println(numpart + xx[ii][0] + " " + xx[ii][1] + " "
						+ xx[numpart][2] + " " + xx[ii][3] + " " + xx[ii][4]
						+ " " + xx[ii][5] + " " + xx[ii][6] + " "
						+ res[numpart]);
			}
		}
		if (ModelName.equals("AdigeCobb3")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				String startDate = "2006-01-01 00:00";
				String endDate = "2007-01-01 00:00";
				int timeStepMinutes = 60;

				OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
				// reader_rain.file =
				// "/Users/giuseppeformetta/Desktop/COBB/totale_pioggeidw_modif2.csv";
				reader_rain.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/KRIGING_RAIN_DISTRIBUTED3bacini.csv";

				// reader_rain.file =
				// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

				reader_rain.idfield = "ID";
				reader_rain.tStart = startDate;
				reader_rain.tEnd = endDate;
				reader_rain.fileNovalue = "-9999.0";
				reader_rain.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
				reader_hydro.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/idrometro.csv";
				reader_hydro.idfield = "Rete_id";
				reader_hydro.tStart = startDate;
				reader_hydro.tEnd = endDate;
				reader_hydro.fileNovalue = "-9999";
				reader_hydro.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
				reader_et.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/ET3bacini.csv";
				reader_et.idfield = "ID";
				reader_et.tStart = startDate;
				reader_et.tEnd = endDate;
				reader_et.fileNovalue = "-9999";
				reader_et.tTimestep = timeStepMinutes;

				int ii = numpart;

				AdigeModified adige = new AdigeModified();
				adige.pCmax = xx[ii][0];
				adige.pB = xx[ii][1];
				adige.pAlpha = xx[ii][2];
				adige.pRs = xx[ii][3];
				adige.pRq = xx[ii][4];
				adige.pMetp = xx[ii][5];
				adige.pMrain = xx[ii][6];
				// adige.pNonlinExp = xx[ii][7];
				adige.pQ0 = 0.2;
				//adige.doNonlinearRes = false;
				adige.pDoReadMeas = true;
				adige.pPathtoMeas = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/MeasuredDischarhe_01012006.csv";
				adige.pDimMeas = 8760;
				adige.pBetanolinear = 1.0;
				adige.pmode = 1;
				adige.pDorouting = false;
				adige.pDoGeom = false;
				adige.fIdDistanze = "dToO";
				adige.fIdnetnum = "netnum";
				adige.inSubbasinDist = basin;
				adige.pm = pm;
				adige.inHillslope = basin;
				adige.fNetnum = "netnum";
				adige.fBaricenter = "avgZ";
				// adige.fVegetation = "uso_reclas";
				adige.inHydrometers = hydrome;
				// adige.inDams = damsFC;
				// adige.inTributary = tributaryFC;
				// adige.inOfftakes = offtakesFC;
				// adige.inVegetation = vegetationData;
				adige.pPfafids = "3";
				// adige.pPfafids = "1";
				adige.fMonpointid = "RETE_ID";
				adige.inNetwork = net;
				adige.fPfaff = "pfafstette";
				adige.fNetelevstart = "elevfirstp";
				adige.fNetelevend = "elevlastpo";

				adige.pRainintensity = -1;
				adige.pRainduration = -1;
				adige.doLog = false;

				adige.tTimestep = timeStepMinutes;
				adige.tStart = startDate;
				adige.tEnd = endDate;
				adige.pNetNumCali = 2;
				reader_rain.initProcess();
				while (reader_rain.doProcess) {

					reader_rain.nextRecord();
					HashMap<Integer, double[]> rainvalues = reader_rain.outData;
					reader_et.nextRecord();
					HashMap<Integer, double[]> etvalues = reader_et.outData;
					adige.inRain = rainvalues;
					adige.inEtp = etvalues;

					reader_hydro.nextRecord();
					adige.inHydrometerdata = reader_hydro.outData;

					adige.process();

				}

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = adige.outMeasured;
				o.simulated = adige.outSimulated;
				o.process();

				// va benissimo
				res[numpart] = 1 - (o.outKGE);
				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				// for(int iii=0;iii<adige.outMeasured.length;iii++){
				// System.out.println( iii+"  "+adige.outMeasured[iii]+ "  "+
				// adige.outSimulated[iii]);
				// }
				System.out.println(numpart + xx[ii][0] + " " + xx[ii][1] + " "
						+ xx[ii][2] + " " + xx[ii][3] + " " + xx[ii][4] + " "
						+ xx[ii][5] + " " + xx[ii][6] + " " + res[numpart]);

			}

		}
		if (ModelName.equals("AdigeCobb1")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				String startDate = "2006-01-01 00:00";
				String endDate = "2007-01-01 00:00";
				int timeStepMinutes = 60;

				OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
				// reader_rain.file =
				// "/Users/giuseppeformetta/Desktop/COBB/totale_pioggeidw_modif2.csv";
				reader_rain.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/KRIGING_RAIN_DISTRIBUTED1bacino.csv";

				// reader_rain.file =
				// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

				reader_rain.idfield = "ID";
				reader_rain.tStart = startDate;
				reader_rain.tEnd = endDate;
				reader_rain.fileNovalue = "-9999.0";
				reader_rain.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
				reader_hydro.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/idrometro.csv";
				reader_hydro.idfield = "Rete_id";
				reader_hydro.tStart = startDate;
				reader_hydro.tEnd = endDate;
				reader_hydro.fileNovalue = "-9999";
				reader_hydro.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
				reader_et.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/ET1bacino.csv";
				reader_et.idfield = "ID";
				reader_et.tStart = startDate;
				reader_et.tEnd = endDate;
				reader_et.fileNovalue = "-9999";
				reader_et.tTimestep = timeStepMinutes;

				int ii = numpart;

				AdigeModified adige = new AdigeModified();
				adige.pCmax = xx[ii][0];
				adige.pB = xx[ii][1];
				adige.pAlpha = xx[ii][2];
				adige.pRs = xx[ii][3];
				adige.pRq = xx[ii][4];
				adige.pMetp = xx[ii][5];
				adige.pMrain = xx[ii][6];
				// adige.pNonlinExp = xx[ii][7];
				adige.pQ0 = 0.0;
				//adige.doNonlinearRes = false;
				adige.pDoReadMeas = true;
				adige.pPathtoMeas = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/MeasuredDischarhe_01012006.csv";
				adige.pDimMeas = 8760;
				adige.pBetanolinear = 1.0;
				adige.pmode = 1;
				adige.pDorouting = false;
				adige.pDoGeom = false;
				adige.fIdDistanze = "dToO";
				adige.fIdnetnum = "netnum";
				adige.inSubbasinDist = basin;
				adige.pm = pm;
				adige.inHillslope = basin;
				adige.fNetnum = "netnum";
				adige.fBaricenter = "avgZ";
				// adige.fVegetation = "uso_reclas";
				adige.inHydrometers = hydrome;
				// adige.inDams = damsFC;
				// adige.inTributary = tributaryFC;
				// adige.inOfftakes = offtakesFC;
				// adige.inVegetation = vegetationData;
				adige.pPfafids = "3";
				// adige.pPfafids = "1";
				adige.fMonpointid = "RETE_ID";
				adige.inNetwork = net;
				adige.fPfaff = "pfafstette";
				adige.fNetelevstart = "elevfirstp";
				adige.fNetelevend = "elevlastpo";

				adige.pRainintensity = -1;
				adige.pRainduration = -1;
				adige.doLog = false;

				adige.tTimestep = timeStepMinutes;
				adige.tStart = startDate;
				adige.tEnd = endDate;
				adige.pNetNumCali = 1;
				reader_rain.initProcess();
				while (reader_rain.doProcess) {

					reader_rain.nextRecord();
					HashMap<Integer, double[]> rainvalues = reader_rain.outData;
					reader_et.nextRecord();
					HashMap<Integer, double[]> etvalues = reader_et.outData;
					adige.inRain = rainvalues;
					adige.inEtp = etvalues;

					reader_hydro.nextRecord();
					adige.inHydrometerdata = reader_hydro.outData;

					adige.process();

				}

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = adige.outMeasured;
				o.simulated = adige.outSimulated;
				o.process();

				// va benissimo
				res[numpart] = 1 - (o.outKGE);
				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				// for(int iii=0;iii<adige.outMeasured.length;iii++){
				// System.out.println( iii+"  "+adige.outMeasured[iii]+ "  "+
				// adige.outSimulated[iii]);
				// }
				System.out.println(numpart + xx[ii][0] + " " + xx[ii][1] + " "
						+ xx[ii][2] + " " + xx[ii][3] + " " + xx[ii][4] + " "
						+ xx[ii][5] + " " + xx[ii][6] + " " + res[numpart]);

			}

		}
		if (ModelName.equals("AdigeCOBBROUTING")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				String startDate = "2006-01-01 00:00";
				String endDate = "2007-01-01 00:00";
				int timeStepMinutes = 60;

				OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
				// reader_rain.file =
				// "/Users/giuseppeformetta/Desktop/COBB/totale_pioggeidw_modif2.csv";
				reader_rain.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/KRIGING_RAIN_DISTRIBUTED.csv";

				// reader_rain.file =
				// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

				reader_rain.idfield = "ID";
				reader_rain.tStart = startDate;
				reader_rain.tEnd = endDate;
				reader_rain.fileNovalue = "-9999.0";
				reader_rain.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
				reader_hydro.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/idrometro.csv";
				reader_hydro.idfield = "Rete_id";
				reader_hydro.tStart = startDate;
				reader_hydro.tEnd = endDate;
				reader_hydro.fileNovalue = "-9999";
				reader_hydro.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
				reader_et.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/totale_ETP_modif.csv";
				reader_et.idfield = "ID";
				reader_et.tStart = startDate;
				reader_et.tEnd = endDate;
				reader_et.fileNovalue = "-9999";
				reader_et.tTimestep = timeStepMinutes;

				int ii = numpart;

				AdigeModified adige = new AdigeModified();
				adige.pCmax = xx[ii][0];
				adige.pB = xx[ii][1];
				adige.pAlpha = xx[ii][2];
				adige.pRs = xx[ii][3];
				adige.pRq = xx[ii][4];
				adige.pMetp = xx[ii][5];
				adige.pMrain = xx[ii][6];
				adige.vR = xx[ii][7];
				;
				adige.pQr = 1.0;
				adige.pLambda1 = xx[ii][8];
				;
				adige.pLambda2 = xx[ii][9];
				;
				adige.pQ0 = 0.1;
				adige.pDoReadMeas = true;
				adige.pPathtoMeas = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/MeasuredDischarhe_01012006.csv";
				adige.pDimMeas = 8760;
				adige.pBetanolinear = 1.0;
				adige.pmode = 1;
				adige.pDorouting = true;
				adige.pDoGeom = false;
				adige.fIdDistanze = "dToO";
				adige.fIdnetnum = "netnum";
				adige.inSubbasinDist = basin;
				adige.pm = pm;
				adige.inHillslope = basin;
				adige.fNetnum = "netnum";
				adige.fBaricenter = "avgZ";
				// adige.fVegetation = "uso_reclas";
				adige.inHydrometers = hydrome;
				// adige.inDams = damsFC;
				// adige.inTributary = tributaryFC;
				// adige.inOfftakes = offtakesFC;
				// adige.inVegetation = vegetationData;
				adige.pPfafids = "3";
				// adige.pPfafids = "1";
				adige.fMonpointid = "RETE_ID";
				adige.inNetwork = net;
				adige.fPfaff = "pfafstette";
				adige.fNetelevstart = "elevfirstp";
				adige.fNetelevend = "elevlastpo";
				adige.pRainintensity = -1;
				adige.pRainduration = -1;
				adige.doLog = false;

				adige.tTimestep = timeStepMinutes;
				adige.tStart = startDate;
				adige.tEnd = endDate;
				adige.pNetNumCali = 4;
				reader_rain.initProcess();
				while (reader_rain.doProcess) {

					reader_rain.nextRecord();
					HashMap<Integer, double[]> rainvalues = reader_rain.outData;
					reader_et.nextRecord();
					HashMap<Integer, double[]> etvalues = reader_et.outData;
					adige.inRain = rainvalues;
					adige.inEtp = etvalues;

					reader_hydro.nextRecord();
					adige.inHydrometerdata = reader_hydro.outData;

					adige.process();

				}

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = adige.outMeasured;
				o.simulated = adige.outSimulated;
				o.process();

				// va benissimo
				res[numpart] = 1 - (o.outKGE);
				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				// for(int iii=0;iii<adige.outMeasured.length;iii++){
				// System.out.println( iii+"  "+adige.outMeasured[iii]+ "  "+
				// adige.outSimulated[iii]);
				// }
				System.out.println(numpart + "  " + xx[ii][0] + " " + xx[ii][1]
						+ " " + xx[ii][2] + " " + xx[ii][3] + " " + xx[ii][4]
						+ " " + xx[ii][5] + " " + xx[ii][6] + " " + xx[ii][7]
						+ " " + xx[ii][8] + " " + xx[ii][9] + " "
						+ res[numpart]);

			}

		}

		/*if (ModelName.equals("AdigeJIMROUTING")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				String startDate = "2005-10-01 00:00";
				String endDate = "2007-10-01 00:00";
				int timeStepMinutes = 60 * 24;

				OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
				reader.file = "/Users/giuseppeformetta/Desktop/JimWORK/Giuseppe/basintempDOK.csv";
				reader.idfield = "ID";
				reader.tStart = startDate;
				reader.tTimestep = 60 * 24;
				reader.tEnd = endDate;
				reader.fileNovalue = "-9999";

				reader.initProcess();

				OmsTimeSeriesIteratorReader reader3 = new OmsTimeSeriesIteratorReader();
				reader3.file = "/Users/giuseppeformetta/Desktop/JimWORK/Jim_dec2012/CondIniI_1_10_2005.csv";
				reader3.idfield = "ID";
				reader3.tStart = startDate;
				reader3.tTimestep = 60 * 24;
				// reader.tEnd = "2000-01-01 00:00";
				reader3.fileNovalue = "-9999";

				reader3.initProcess();

				OmsTimeSeriesIteratorReader reader4 = new OmsTimeSeriesIteratorReader();
				reader4.file = "/Users/giuseppeformetta/Desktop/JimWORK/Jim_dec2012/CondIniL_1_10_2005.csv";
				reader4.idfield = "ID";
				reader4.tStart = startDate;
				reader4.tTimestep = 60 * 24;
				// reader.tEnd = "2000-01-01 00:00";
				reader4.fileNovalue = "-9999";

				reader4.initProcess();

				OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
				// reader_rain.file =
				// "/Users/giuseppeformetta/Desktop/COBB/totale_pioggeidw_modif2.csv";
				reader_rain.file = "/Users/giuseppeformetta/Desktop/JimWORK/Giuseppe/basinRainfallDOK.csv";

				// reader_rain.file =
				// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

				reader_rain.idfield = "ID";
				reader_rain.tStart = startDate;
				reader_rain.tEnd = endDate;
				reader_rain.fileNovalue = "-9999.0";
				reader_rain.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
				reader_hydro.file = "/Users/giuseppeformetta/Desktop/JimWORK/Jim_dec2012/hydro.csv";
				reader_hydro.idfield = "ID";
				reader_hydro.tStart = startDate;
				reader_hydro.tEnd = endDate;
				reader_hydro.fileNovalue = "-9999";
				reader_hydro.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
				reader_et.file = "/Users/giuseppeformetta/Desktop/JimWORK/Giuseppe/BasinETP.csv";
				reader_et.idfield = "ID";
				reader_et.tStart = startDate;
				reader_et.tEnd = endDate;
				reader_et.fileNovalue = "-9999";
				reader_et.tTimestep = timeStepMinutes;

				int ii = numpart;

				SnowMeltingUpdated28012013 insolation = new SnowMeltingUpdated28012013();
				insolation.inDem = inDem;
				insolation.inStations = inStations;
				insolation.fStationsid = "netnum";
				insolation.inInsFeb = inInsFeb;
				insolation.inInsMar = inInsMar;
				insolation.inInsApr = inInsApr;
				insolation.inInsMay = inInsMay;
				insolation.inInsJun = inInsJun;
				insolation.inSkyview = inSkyview;
				// insolation.defaultLapse=-.0065;
				// insolation.defaultRH=0.4;
				// insolation.defaultVisibility=60;
				insolation.doRaster = false;
				insolation.pCmf = xx[numpart][7];
				insolation.pCff = xx[numpart][8];
				insolation.pR = xx[numpart][9];
				insolation.pCr = xx[numpart][10];
				insolation.pCs = xx[numpart][11];
				insolation.pTmelt = xx[numpart][12];
				insolation.doDaily = true;

				AdigeModified adige = new AdigeModified();
				adige.pCmax = xx[ii][0];
				adige.pB = xx[ii][1];
				adige.pAlpha = xx[ii][2];
				adige.pRs = xx[ii][3];
				adige.pRq = xx[ii][4];
				adige.pMetp = xx[ii][6];
				adige.pMrain = xx[ii][5];
				// adige.vR = xx[ii][7];
				;
				adige.pQr = 1.0;
				// adige.pLambda1 = xx[ii][8];
				;
				// adige.pLambda2 = xx[ii][9];
				;
				adige.pQ0 = 0.0;
				adige.pDoReadMeas = true;
				adige.pPathtoMeas = "/Users/giuseppeformetta/Desktop/JimWORK/DischarheJimUSGSOKKKKK2005.txt";
				adige.pDimMeas = 365*2;
				adige.pBetanolinear = 1.0;
				adige.pmode = 1;
				adige.pDorouting = false;
				adige.pDoGeom = false;
				adige.fIdDistanze = "dToO";
				adige.fIdnetnum = "netnum";
				adige.inSubbasinDist = basin;
				adige.pm = pm;
				adige.inHillslope = basin;
				adige.fNetnum = "netnum";
				adige.fBaricenter = "avgZ";
				// adige.fVegetation = "uso_reclas";
				adige.inHydrometers = hydrome;
				// adige.inDams = damsFC;
				// adige.inTributary = tributaryFC;
				// adige.inOfftakes = offtakesFC;
				// adige.inVegetation = vegetationData;
				adige.pPfafids = "1";
				// adige.pPfafids = "1";
				adige.fMonpointid = "SITE_NUMBE";
				adige.inNetwork = net;
				adige.fPfaff = "pfafstette";
				adige.fNetelevstart = "elevfirstp";
				adige.fNetelevend = "elevlastpo";
				adige.pRainintensity = -1;
				adige.pRainduration = -1;
				adige.doLog = false;

				adige.tTimestep = timeStepMinutes;
				adige.tStart = startDate;
				adige.tEnd = endDate;
				adige.pNetNumCali = 6;
				reader_rain.initProcess();
				boolean flag = true;

				while (reader_rain.doProcess) {
					reader.nextRecord();
					HashMap<Integer, double[]> temp = reader.outData;
					reader_rain.nextRecord();
					HashMap<Integer, double[]> rainvalues = reader_rain.outData;
					reader_et.nextRecord();
					HashMap<Integer, double[]> etvalues = reader_et.outData;
					// ///////////////////NEVE////////////////////////////////////////
					insolation.time = reader.tCurrent;

					insolation.inTemp = temp;
					insolation.inRainfall = rainvalues;
					reader3.nextRecord();
					reader4.nextRecord();

					if (flag == true) {
						insolation.inCondIniI = reader3.outData;
						insolation.inCondIniL = reader4.outData;
						flag = false;
					}
					insolation.process();
					// ///////////////////NEVE////////////////////////////////////////

					adige.inRain = rainvalues;
					adige.inEtp = etvalues;
					adige.inMelting = insolation.outMeltingData;
					reader_hydro.nextRecord();
					adige.inHydrometerdata = reader_hydro.outData;

					adige.process();

				}

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = adige.outMeasured;
				o.simulated = adige.outSimulated;
				o.process();

				// va benissimo
				res[numpart] = 1 - (o.outKGE);
				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				// for(int iii=0;iii<adige.outMeasured.length;iii++){
				// System.out.println( iii+"  "+adige.outMeasured[iii]+
				// "  "+
				// adige.outSimulated[iii]);
				// }
				System.out.println(numpart + "  " + xx[ii][0] + " " + xx[ii][1]
						+ " " + xx[ii][2] + " " + xx[ii][3] + " " + xx[ii][4]
						+ " " + xx[ii][5] + " " + xx[ii][6] + " " + xx[ii][7]
						+ " " + xx[ii][8] + " " + xx[ii][9] + " " + xx[ii][10]
						+ " " + xx[ii][11] + " " + +res[numpart]);

			}

		}*/
		/*if (ModelName.equals("AdigeJIMNOROUTING")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				String startDate = "2005-05-13 00:00";
				String endDate = "2008-05-13 00:00";
				int timeStepMinutes = 60 * 24;

				OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
				reader.file = "/Users/giuseppeformetta/Desktop/JimWORK/Jim_dec2012/JAMI_JIM/BasinTmean.csv";
				reader.idfield = "ID";
				reader.tStart = startDate;
				reader.tTimestep = 60 * 24;
				reader.tEnd = endDate;
				reader.fileNovalue = "-9999";

				reader.initProcess();

				OmsTimeSeriesIteratorReader reader3 = new OmsTimeSeriesIteratorReader();
				reader3.file = "/Users/giuseppeformetta/Desktop/JimWORK/Jim_dec2012/CondIniI.csv";
				reader3.idfield = "ID";
				reader3.tStart = startDate;
				reader3.tTimestep = 60 * 24;
				// reader.tEnd = "2000-01-01 00:00";
				reader3.fileNovalue = "-9999";

				reader3.initProcess();

				OmsTimeSeriesIteratorReader reader4 = new OmsTimeSeriesIteratorReader();
				reader4.file = "/Users/giuseppeformetta/Desktop/JimWORK/Jim_dec2012/CondIniL.csv";
				reader4.idfield = "ID";
				reader4.tStart = startDate;
				reader4.tTimestep = 60 * 24;
				// reader.tEnd = "2000-01-01 00:00";
				reader4.fileNovalue = "-9999";

				reader4.initProcess();

				OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
				// reader_rain.file =
				// "/Users/giuseppeformetta/Desktop/COBB/totale_pioggeidw_modif2.csv";
				reader_rain.file = "/Users/giuseppeformetta/Desktop/JimWORK/Jim_dec2012/basinRainfall.csv";

				// reader_rain.file =
				// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

				reader_rain.idfield = "ID";
				reader_rain.tStart = startDate;
				reader_rain.tEnd = endDate;
				reader_rain.fileNovalue = "-9999.0";
				reader_rain.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
				reader_hydro.file = "/Users/giuseppeformetta/Desktop/JimWORK/Jim_dec2012/hydro.csv";
				reader_hydro.idfield = "ID";
				reader_hydro.tStart = startDate;
				reader_hydro.tEnd = endDate;
				reader_hydro.fileNovalue = "-9999";
				reader_hydro.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
				reader_et.file = "/Users/giuseppeformetta/Desktop/JimWORK/Jim_dec2012/BasinETP.csv";
				reader_et.idfield = "ID";
				reader_et.tStart = startDate;
				reader_et.tEnd = endDate;
				reader_et.fileNovalue = "-9999";
				reader_et.tTimestep = timeStepMinutes;

				int ii = numpart;

				SnowMeltingUpdated28012013 insolation = new SnowMeltingUpdated28012013();
				insolation.inDem = inDem;
				insolation.inStations = inStations;
				insolation.fStationsid = "netnum";
				insolation.inInsFeb = inInsFeb;
				insolation.inInsMar = inInsMar;
				insolation.inInsApr = inInsApr;
				insolation.inInsMay = inInsMay;
				insolation.inInsJun = inInsJun;
				insolation.inSkyview = inSkyview;
				// insolation.defaultLapse=-.0065;
				// insolation.defaultRH=0.4;
				// insolation.defaultVisibility=60;
				insolation.doRaster = false;
				insolation.pCmf = xx[numpart][10];
				insolation.pCff = xx[numpart][11];
				insolation.pR = xx[numpart][12];
				insolation.pCr = 1.0;
				insolation.pCs = 1.0;
				insolation.pTmelt = 1.0;
				insolation.doDaily = true;

				AdigeModified adige = new AdigeModified();
				adige.pCmax = xx[ii][0];
				adige.pB = xx[ii][1];
				adige.pAlpha = xx[ii][2];
				adige.pRs = xx[ii][3];
				adige.pRq = xx[ii][4];
				adige.pMetp = xx[ii][5];
				adige.pMrain = xx[ii][6];
				adige.vR = xx[ii][7];
				adige.pQr = 1.0;
				adige.pLambda1 = xx[ii][8];
				adige.pLambda2 = xx[ii][9];

				adige.pQ0 = 0.1;
				adige.pDoReadMeas = true;
				adige.pPathtoMeas = "/Users/giuseppeformetta/Desktop/JimWORK/Portate.txt";
				adige.pDimMeas = 977;
				adige.pBetanolinear = 1.0;
				adige.pmode = 1;
				adige.pDorouting = true;
				adige.pDoGeom = false;
				adige.fIdDistanze = "dToO";
				adige.fIdnetnum = "netnum";
				adige.inSubbasinDist = basin;
				adige.pm = pm;
				adige.inHillslope = basin;
				adige.fNetnum = "netnum";
				adige.fBaricenter = "avgZ";
				// adige.fVegetation = "uso_reclas";
				adige.inHydrometers = hydrome;
				// adige.inDams = damsFC;
				// adige.inTributary = tributaryFC;
				// adige.inOfftakes = offtakesFC;
				// adige.inVegetation = vegetationData;
				adige.pPfafids = "3";
				// adige.pPfafids = "1";
				adige.fMonpointid = "POINTID";
				adige.inNetwork = net;
				adige.fPfaff = "pfafstette";
				adige.fNetelevstart = "elevfirstp";
				adige.fNetelevend = "elevlastpo";
				adige.pRainintensity = -1;
				adige.pRainduration = -1;
				adige.doLog = false;

				adige.tTimestep = timeStepMinutes;
				adige.tStart = startDate;
				adige.tEnd = endDate;
				adige.pNetNumCali = 5;
				reader_rain.initProcess();
				boolean flag = true;

				while (reader_rain.doProcess) {
					reader.nextRecord();
					HashMap<Integer, double[]> temp = reader.outData;
					reader_rain.nextRecord();
					HashMap<Integer, double[]> rainvalues = reader_rain.outData;
					reader_et.nextRecord();
					HashMap<Integer, double[]> etvalues = reader_et.outData;
					// ///////////////////NEVE////////////////////////////////////////
					insolation.time = reader.tCurrent;

					HashMap<Integer, double[]> pio = reader_rain.outData;
					insolation.inTemp = temp;
					insolation.inRainfall = pio;
					reader3.nextRecord();
					reader4.nextRecord();

					if (flag == true) {
						insolation.inCondIniI = reader3.outData;
						insolation.inCondIniL = reader4.outData;
						flag = false;
					}
					insolation.process();
					// ///////////////////NEVE////////////////////////////////////////

					adige.inRain = rainvalues;
					adige.inEtp = etvalues;
					adige.inMelting = insolation.outMeltingData;
					reader_hydro.nextRecord();
					adige.inHydrometerdata = reader_hydro.outData;

					adige.process();

				}

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = adige.outMeasured;
				o.simulated = adige.outSimulated;
				o.process();

				// va benissimo
				res[numpart] = 1 - (o.outKGE);
				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				// for(int iii=0;iii<adige.outMeasured.length;iii++){
				// System.out.println( iii+"  "+adige.outMeasured[iii]+
				// "  "+
				// adige.outSimulated[iii]);
				// }
				System.out.println(numpart + "  " + xx[ii][0] + " " + xx[ii][1]
						+ " " + xx[ii][2] + " " + xx[ii][3] + " " + xx[ii][4]
						+ " " + xx[ii][5] + " " + xx[ii][6] + " " + xx[ii][7]
						+ " " + xx[ii][8] + " " + xx[ii][9] + " "
						+ res[numpart]);

			}

		}
		if (ModelName.equals("AdigeCOBBROUTING3")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				String startDate = "2006-01-01 00:00";
				String endDate = "2007-01-01 00:00";
				int timeStepMinutes = 60;

				TimeSeriesIteratorReader reader_rain = new TimeSeriesIteratorReader();
				// reader_rain.file =
				// "/Users/giuseppeformetta/Desktop/COBB/totale_pioggeidw_modif2.csv";
				reader_rain.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/Piogge3bacini.csv";

				// reader_rain.file =
				// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

				reader_rain.idfield = "ID";
				reader_rain.tStart = startDate;
				reader_rain.tEnd = endDate;
				reader_rain.fileNovalue = "-9999.0";
				reader_rain.tTimestep = timeStepMinutes;

				TimeSeriesIteratorReader reader_hydro = new TimeSeriesIteratorReader();
				reader_hydro.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/idrometro.csv";
				reader_hydro.idfield = "Rete_id";
				reader_hydro.tStart = startDate;
				reader_hydro.tEnd = endDate;
				reader_hydro.fileNovalue = "-9999";
				reader_hydro.tTimestep = timeStepMinutes;

				TimeSeriesIteratorReader reader_et = new TimeSeriesIteratorReader();
				reader_et.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/ET3bacini.csv";
				reader_et.idfield = "ID";
				reader_et.tStart = startDate;
				reader_et.tEnd = endDate;
				reader_et.fileNovalue = "-9999";
				reader_et.tTimestep = timeStepMinutes;

				int ii = numpart;

				AdigeModified adige = new AdigeModified();
				adige.pCmax = xx[ii][0];
				adige.pB = xx[ii][1];
				adige.pAlpha = xx[ii][2];
				adige.pRs = xx[ii][3];
				adige.pRq = xx[ii][4];
				adige.pMetp = xx[ii][5];
				adige.pMrain = xx[ii][6];
				adige.vR = xx[ii][7];
				;
				adige.pQr = 1.0;
				adige.pLambda1 = xx[ii][8];
				;
				adige.pLambda2 = xx[ii][9];
				;
				adige.pQ0 = 0.1;
				adige.pDoReadMeas = true;
				adige.pPathtoMeas = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/MeasuredDischarhe_01012006.csv";
				adige.pDimMeas = 8760;
				adige.pBetanolinear = 1.0;
				adige.pmode = 1;
				adige.pDorouting = true;
				adige.pDoGeom = false;
				adige.fIdDistanze = "dToO";
				adige.fIdnetnum = "netnum";
				adige.inSubbasinDist = basin;
				adige.pm = pm;
				adige.inHillslope = basin;
				adige.fNetnum = "netnum";
				adige.fBaricenter = "avgZ";
				// adige.fVegetation = "uso_reclas";
				adige.inHydrometers = hydrome;
				// adige.inDams = damsFC;
				// adige.inTributary = tributaryFC;
				// adige.inOfftakes = offtakesFC;
				// adige.inVegetation = vegetationData;
				adige.pPfafids = "3";
				// adige.pPfafids = "1";
				adige.fMonpointid = "RETE_ID";
				adige.inNetwork = net;
				adige.fPfaff = "pfafstette";
				adige.fNetelevstart = "elevfirstp";
				adige.fNetelevend = "elevlastpo";
				adige.pRainintensity = -1;
				adige.pRainduration = -1;
				adige.doLog = false;

				adige.tTimestep = timeStepMinutes;
				adige.tStart = startDate;
				adige.tEnd = endDate;
				adige.pNetNumCali = 2;
				reader_rain.initProcess();
				while (reader_rain.doProcess) {

					reader_rain.nextRecord();
					HashMap<Integer, double[]> rainvalues = reader_rain.outData;
					reader_et.nextRecord();
					HashMap<Integer, double[]> etvalues = reader_et.outData;
					adige.inRain = rainvalues;
					adige.inEtp = etvalues;

					reader_hydro.nextRecord();
					adige.inHydrometerdata = reader_hydro.outData;

					adige.process();

				}

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = adige.outMeasured;
				o.simulated = adige.outSimulated;
				o.process();

				// va benissimo
				res[numpart] = 1 - (o.outKGE);
				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				// for(int iii=0;iii<adige.outMeasured.length;iii++){
				// System.out.println( iii+"  "+adige.outMeasured[iii]+ "  "+
				// adige.outSimulated[iii]);
				// }
				System.out.println(numpart + "  " + xx[ii][0] + " " + xx[ii][1]
						+ " " + xx[ii][2] + " " + xx[ii][3] + " " + xx[ii][4]
						+ " " + xx[ii][5] + " " + xx[ii][6] + " " + xx[ii][7]
						+ " " + xx[ii][8] + " " + xx[ii][9] + " "
						+ res[numpart]);

			}

		}*/
		if (ModelName.equals("AdigeCOBBROUTING1")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				String startDate = "2006-01-01 00:00";
				String endDate = "2007-01-01 00:00";
				int timeStepMinutes = 60;

				OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
				// reader_rain.file =
				// "/Users/giuseppeformetta/Desktop/COBB/totale_pioggeidw_modif2.csv";
				reader_rain.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/KRIGING_RAIN_DISTRIBUTED1bacino.csv";

				// reader_rain.file =
				// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

				reader_rain.idfield = "ID";
				reader_rain.tStart = startDate;
				reader_rain.tEnd = endDate;
				reader_rain.fileNovalue = "-9999.0";
				reader_rain.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
				reader_hydro.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/idrometro.csv";
				reader_hydro.idfield = "Rete_id";
				reader_hydro.tStart = startDate;
				reader_hydro.tEnd = endDate;
				reader_hydro.fileNovalue = "-9999";
				reader_hydro.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
				reader_et.file = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/ET1bacino.csv";
				reader_et.idfield = "ID";
				reader_et.tStart = startDate;
				reader_et.tEnd = endDate;
				reader_et.fileNovalue = "-9999";
				reader_et.tTimestep = timeStepMinutes;

				int ii = numpart;

				AdigeModified adige = new AdigeModified();
				adige.pCmax = xx[ii][0];
				adige.pB = xx[ii][1];
				adige.pAlpha = xx[ii][2];
				adige.pRs = xx[ii][3];
				adige.pRq = xx[ii][4];
				adige.pMetp = xx[ii][5];
				adige.pMrain = xx[ii][6];
				adige.vR = xx[ii][7];
				;
				adige.pQr = 1.0;
				adige.pLambda1 = xx[ii][8];
				;
				adige.pLambda2 = xx[ii][9];
				;
				adige.pQ0 = 0.1;
				adige.pDoReadMeas = true;
				adige.pPathtoMeas = "/Users/giuseppeformetta/Desktop/paperEMS/COBB/MeasuredDischarhe_01012006.csv";
				adige.pDimMeas = 8760;
				adige.pBetanolinear = 1.0;
				adige.pmode = 1;
				adige.pDorouting = true;
				adige.pDoGeom = false;
				adige.fIdDistanze = "dToO";
				adige.fIdnetnum = "netnum";
				adige.inSubbasinDist = basin;
				adige.pm = pm;
				adige.inHillslope = basin;
				adige.fNetnum = "netnum";
				adige.fBaricenter = "avgZ";
				// adige.fVegetation = "uso_reclas";
				adige.inHydrometers = hydrome;
				// adige.inDams = damsFC;
				// adige.inTributary = tributaryFC;
				// adige.inOfftakes = offtakesFC;
				// adige.inVegetation = vegetationData;
				adige.pPfafids = "3";
				// adige.pPfafids = "1";
				adige.fMonpointid = "RETE_ID";
				adige.inNetwork = net;
				adige.fPfaff = "pfafstette";
				adige.fNetelevstart = "elevfirstp";
				adige.fNetelevend = "elevlastpo";
				adige.pRainintensity = -1;
				adige.pRainduration = -1;
				adige.doLog = false;

				adige.tTimestep = timeStepMinutes;
				adige.tStart = startDate;
				adige.tEnd = endDate;
				adige.pNetNumCali = 1;
				reader_rain.initProcess();
				while (reader_rain.doProcess) {

					reader_rain.nextRecord();
					HashMap<Integer, double[]> rainvalues = reader_rain.outData;
					reader_et.nextRecord();
					HashMap<Integer, double[]> etvalues = reader_et.outData;
					adige.inRain = rainvalues;
					adige.inEtp = etvalues;

					reader_hydro.nextRecord();
					adige.inHydrometerdata = reader_hydro.outData;

					adige.process();

				}

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = adige.outMeasured;
				o.simulated = adige.outSimulated;
				o.process();

				// va benissimo
				res[numpart] = 1 - (o.outKGE);
				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				// for(int iii=0;iii<adige.outMeasured.length;iii++){
				// System.out.println( iii+"  "+adige.outMeasured[iii]+ "  "+
				// adige.outSimulated[iii]);
				// }
				System.out.println(numpart + "  " + xx[ii][0] + " " + xx[ii][1]
						+ " " + xx[ii][2] + " " + xx[ii][3] + " " + xx[ii][4]
						+ " " + xx[ii][5] + " " + xx[ii][6] + " " + xx[ii][7]
						+ " " + xx[ii][8] + " " + xx[ii][9] + " "
						+ res[numpart]);

			}

		}
		if (ModelName.equals("AdigeVeneto")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				String startDate = "2002-10-01 00:00";
				String endDate = "2003-07-21 15:00";
				int timeStepMinutes = 60;

				OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
				reader_rain.file = "/Users/giuseppeformetta/Desktop/VENETO_MATTEO/csv/PIO_ok.csv";
				// reader_rain.file =
				// "/Users/formeppe/Desktop/ORARIO/pioggeOK_orarie_2002_2008_NEWFINELESS.csv";

				reader_rain.idfield = "ID";
				reader_rain.tStart = startDate;
				reader_rain.tEnd = endDate;
				reader_rain.fileNovalue = "-9999.0";
				reader_rain.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_hydro = new OmsTimeSeriesIteratorReader();
				reader_hydro.file = "/Users/giuseppeformetta/Desktop/VENETO_MATTEO/csv/hydro.csv";
				reader_hydro.idfield = "ID";
				reader_hydro.tStart = startDate;
				reader_hydro.tEnd = endDate;
				reader_hydro.fileNovalue = "-9999";
				reader_hydro.tTimestep = timeStepMinutes;

				OmsTimeSeriesIteratorReader reader_et = new OmsTimeSeriesIteratorReader();
				reader_et.file = "/Users/giuseppeformetta/Desktop/VENETO_MATTEO/csv/ET_ok.csv";
				reader_et.idfield = "ID";
				reader_et.tStart = startDate;
				reader_et.tEnd = endDate;
				reader_et.fileNovalue = "-9999";
				reader_et.tTimestep = timeStepMinutes;

				int ii = numpart;

				AdigeModified adige = new AdigeModified();
				adige.pCmax = xx[ii][0];
				adige.pB = xx[ii][1];
				adige.pAlpha = xx[ii][2];
				adige.pRs = xx[ii][3];
				adige.pRq = xx[ii][4];
				adige.pMetp = xx[ii][5];
				adige.pMrain = xx[ii][6];
				adige.pQ0 = 0.0;
				adige.pDoReadMeas = true;
				adige.pPathtoMeas = "/Users/giuseppeformetta/Desktop/VENETO_MATTEO/csv/dok";
				adige.pDimMeas = 9960;
				adige.pBetanolinear = 1.0;
				adige.pmode = 1;
				adige.pDorouting = false;
				adige.pDoGeom = true;
				adige.fIdDistanze = "dist";
				adige.fIdnetnum = "netnum";
				adige.inSubbasinDist = basin;
				adige.pm = pm;
				adige.inHillslope = basin;
				adige.fNetnum = "netnum";
				adige.fBaricenter = "avgZ";
				// adige.fVegetation = "uso_reclas";
				adige.inHydrometers = hydrome;
				// adige.inDams = damsFC;
				// adige.inTributary = tributaryFC;
				// adige.inOfftakes = offtakesFC;
				// adige.inVegetation = vegetationData;
				adige.pPfafids = "3";
				// adige.pPfafids = "1";
				adige.fMonpointid = "ID";
				adige.inNetwork = net;
				adige.fPfaff = "pfafstette";
				adige.fNetelevstart = "elevfirstp";
				adige.fNetelevend = "elevlastpo";

				adige.pRainintensity = -1;
				adige.pRainduration = -1;
				adige.doLog = false;

				adige.tTimestep = timeStepMinutes;
				adige.tStart = startDate;
				adige.tEnd = endDate;
				adige.pNetNumCali = 4;
				reader_rain.initProcess();
				while (reader_rain.doProcess) {

					reader_rain.nextRecord();
					HashMap<Integer, double[]> rainvalues = reader_rain.outData;
					reader_et.nextRecord();
					HashMap<Integer, double[]> etvalues = reader_et.outData;
					adige.inRain = rainvalues;
					adige.inEtp = etvalues;

					reader_hydro.nextRecord();
					adige.inHydrometerdata = reader_hydro.outData;

					adige.process();

				}

				reader_rain.close();
				reader_et.close();
				reader_hydro.close();

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = adige.outMeasured;
				o.simulated = adige.outSimulated;
				o.process();
				// va benissimo
				res[numpart] = 1 - (o.outKGE);
				// res[numpart] = Math.pow(((1 - o.outFHF) * (o.outFHF) +
				// o.outFLF
				// * o.outFLF), 0.5);
				// res[numpart] = o.outObFunRMSE;
				System.out.println(numpart + "    " + res[numpart]);
				System.out.println(numpart + "    " + res[numpart]);

			}

		}
		if (ModelName.equals("Vgm")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {

				VGM v = new VGM();
				v.modelname = inmodelname;
				v.distances = inDistance;
				v.nugget = xx[numpart][0];
				v.sill = xx[numpart][1];
				v.Np = inNp;
				v.range = xx[numpart][2];
				v.process();

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = inVar;
				o.simulated = v.result;
				o.w = v.Npairs;
				o.process();
				// va benissimo
				// res[numpart] =
				// Math.pow((o.outObFunRMSEW)*(o.outObFunRMSEW)+(1-o.outKGE)*(1-o.outKGE),0.5);

				res[numpart] = (o.outObFunRMSE);
				// res[numpart] = (o.outKGE);

			}

		}
		if (ModelName.equals("PeakFlow")) {

			for (int numpart = 0; numpart < xx.length; numpart++) {

				OmsRescaledDistance rescaledDistance_10 = new OmsRescaledDistance();
				rescaledDistance_10.inFlow = inDrainDir;
				rescaledDistance_10.inNet = inNet;
				rescaledDistance_10.pRatio = xx[numpart][0];
				rescaledDistance_10.pm = pm;

				rescaledDistance_10.process();

				GridCoverage2D rescaledDistance10 = rescaledDistance_10.outRescaled;

				OmsRescaledDistance rescaledDistance_100 = new OmsRescaledDistance();
				rescaledDistance_100.inFlow = inDrainDir;
				rescaledDistance_100.inNet = inNet;
				rescaledDistance_100.pRatio = xx[numpart][1];
				;
				rescaledDistance_100.pm = pm;

				rescaledDistance_100.process();

				GridCoverage2D rescaledDistance100 = rescaledDistance_100.outRescaled;

				//

				OmsTimeSeriesReader reader = new OmsTimeSeriesReader();
				reader.file = pathToRain;
				reader.read();
				HashMap<DateTime, double[]> outData = reader.outData;

				OmsPeakflow pp = new OmsPeakflow();
				pp.inRescaledsub = rescaledDistance100;
				pp.inRescaledsup = rescaledDistance10;
				pp.inTopindex = inTopIndex;
				pp.pCelerity = xx[numpart][2];

				pp.pDiffusion = xx[numpart][3];

				pp.inRainfall = outData;
				pp.pSat = xx[numpart][4];

			//	pp.pMrain = xx[numpart][5];

				pp.outputStepArg = 3600;
				pp.process();
				HashMap<DateTime, double[]> a = pp.outDischarge;

				OmsTimeSeriesReader reader2 = new OmsTimeSeriesReader();
				reader2.file = pathToRunoff;
				reader2.read();
				HashMap<DateTime, double[]> outData2 = reader2.outData;
				double[] misurati = new double[outData2.size()];
				double[] modellati = new double[outData2.size()];
				int cout = 0;
				Set<Entry<DateTime, double[]>> entrySet = outData2.entrySet();
				Set<Entry<DateTime, double[]>> entrySet2 = pp.outDischarge
						.entrySet();

				for (Entry<DateTime, double[]> entry : entrySet) {
					DateTime dateTime = entry.getKey();
					double[] valuesArray = entry.getValue();
					misurati[cout] = valuesArray[0];
					cout += 1;
				}
				cout = 0;
				for (Entry<DateTime, double[]> entry : entrySet2) {
					DateTime dateTime = entry.getKey();
					double[] valuesArray = entry.getValue();
					modellati[cout] = valuesArray[0];
					cout += 1;
					System.out.println(cout);
					if (cout >= misurati.length) {
						break;
					}
				}

				ObFun o = new ObFun();
				o.observed = misurati;
				o.simulated = modellati;
				o.process();
				res[numpart] = 1.0 - (o.outKGE);
				System.out.println("KGE= " + o.outKGE + "  NSE= "
						+ o.outObFunNash);
			}
		}
		/*if (ModelName.equals("HyMod")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				HyMod prova = new HyMod();

				prova.Prec = new Matrix(Extra_Precip, int_Extra_MaxT);
				prova.ETP = new Matrix(Extra_PET, int_Extra_MaxT);
				prova.Cmax = xx[numpart][0];
				prova.b = xx[numpart][1];
				prova.alpha = xx[numpart][2];
				prova.Rs = xx[numpart][3];
				prova.Rq = xx[numpart][4];
				prova.Area = 11;
				prova.timestep = 24.0;
				prova.Q0 = xx[numpart][5];

				prova.process();
				Matrix ModelPred_m = new Matrix(int_Extra_MaxT - Warmup, 1);
				ModelPred_m = prova.outHyMod;
				double[][] ModelPred = new double[int_Extra_MaxT - Warmup][1];
				ModelPred = ModelPred_m.getArray();
				// tutti i dati
				// double[] modpredvect=new double [ModelPred.length];
				// double[] observedvect=new double [ModelPred.length];
				// dati vrugt
				double[] modpredvect = new double[ModelPred.length - Warmup + 1];
				double[] observedvect = new double[ModelPred.length - Warmup
						+ 1];

				// int conta = -1;
				// for (int kk = Warmup; kk < ModelPred.length; kk++) {
				// if (Measurement_MeasData[kk - Warmup] != -999
				// && ModelPred[kk - Warmup][0] != -999) {
				// // if(Measurement_MeasData[kk]!=-999 &&
				// // ModelPred[kk][0]!=-999){
				//
				// conta += 1;
				// // per tutti i dati
				// // observedvect[conta]=Measurement_MeasData[kk];
				// // modpredvect[conta]= ModelPred[kk-Warmup][0];
				// // modpredvect[conta]= ModelPred[kk][0];
				// // per dream di vrugt
				// observedvect[conta] = Measurement_MeasData[kk - Warmup];
				// modpredvect[conta] = ModelPred[kk - Warmup][0];
				// System.out.println(conta+"    "+observedvect[conta]+"    "+modpredvect[conta]);
				// }
				// }
				for (int kk = Warmup; kk < ModelPred.length - Warmup; kk++) {
					// System.out.println(kk);
					observedvect[kk] = Measurement_MeasData[kk + Warmup - 1];
					modpredvect[kk] = ModelPred[kk][0];
				}

				ObFun o = new ObFun();
				o.NaN = -999;
				o.observed = observedvect;
				o.simulated = modpredvect;
				o.process();
				// va benissimo
				res[numpart] = 1.0 - o.outKGE;

				// res[numpart] = 1 - (o.outKGE);
				// res[numpart] = 1 - o.outObFunNash;
				// res[numpart]=Math.pow(((1-o.outFHF)*(o.outFHF)+o.outFLF*o.outFLF),0.5);

				// matrice[dimmat][0]=o.outObFunNash;
				// matrice[dimmat][1]=o.outKGE;
				// dimmat++;
				// System.out.println(dimmat+"   NSE="+o.outObFunNash+" FLS="+o.outObFunFsls+" PBIAS="+o.outObFunPBIAS+" RMSE="+o.outObFunRMSE);
				// System.out.println("Cmax= "+ xx[numpart][0]);
				// System.out.println("b= "+ xx[numpart][1]);
				// System.out.println("a= "+ xx[numpart][2]);
				// System.out.println("rs= "+ xx[numpart][3]);
				// System.out.println("rq= "+ xx[numpart][4]);
				//
				//
				System.out.println("NSE=" + o.outObFunNash + " FLS="
						+ o.outObFunFsls + " PBIAS=" + o.outObFunPBIAS
						+ " RMSE=" + o.outObFunRMSE);
			}
		}*/
		
		
		return res;

	}

	public double[] Compute_gBest(double xx[][], double[] vettcostnew) {
		double re[] = g_best;
		int pos = 0;
		double min = g_best_value;
		for (int i = 0; i < vettcostnew.length; i++) {
			if (Math.abs(vettcostnew[i]) <= min) {
				g_best_value = Math.abs(vettcostnew[i]);
				min = Math.abs(vettcostnew[i]);
				pos = i;
				for (int ii = 0; ii < xx[0].length; ii++) {
					re[ii] = xx[pos][ii];
				}
			}
		}

		// System.out.println("minimo="+g_best_value);
		return re;
	}

	public double[][] Compute_particle(double pos[][], double[][] vel) {

		double xnew[][] = new double[pos.length][pos[0].length];
		for (int i = 0; i < vel.length; i++) {
			for (int j = 0; j < vel[0].length; j++) {
				xnew[i][j] = pos[i][j] + vel[i][j];
			}

		}
		double[][] xneww = ReflectBounds(xnew, ParRange_maxn, ParRange_minn);
		return xneww;
	}

	public double[][] Compute_velocity(double pos[][], double[][] vel) {

		double velnew[][] = new double[pos.length][pos[0].length];
		for (int i = 0; i < vel.length; i++) {
			for (int j = 0; j < vel[0].length; j++) {
				double c1 = 1.5;
				double r1 = Math.random();
				double c2 = 2.5;
				double r2 = Math.random();
				double inertia = 0.5;
				velnew[i][j] = inertia * vel[i][j] + c1 * r1
						* (p_best[i][j] - pos[i][j]) + c2 * r2
						* (g_best[j] - pos[i][j]);
			}

		}
		return velnew;
	}

	public double[][] Compute_pBest(double currentpos[][], double[] currentbest) {
		double pos_best[][] = p_best;
		for (int i = 0; i < currentbest.length; i++) {
			// per tutti
			if (Math.abs(currentbest[i]) < Math.abs(p_best_value[i])) {
				// per nash
				// if(Math.abs(currentbest[i])>Math.abs(p_best_value[i])){
				p_best_value[i] = Math.abs(currentbest[i]);
				for (int j = 0; j < currentpos[0].length; j++) {
					pos_best[i][j] = currentpos[i][j];
				}
			}
		}
		return pos_best;
	}

	public double[][] ReflectBounds(double[][] neww, double[] ParRange_maxnn,
			double[] ParRange_minnn) {
		// ParRange_maxnn=new double []{15000,50,0.99,0.3,0.9};
		// ParRange_minnn=new double []{1,0.01,0.01,0.0001,0.01};
		// Checks the bounds of the parameters
		// First determine the size of new
		// int nmbOfIndivs=neww.length;
		// int Dim=neww[0].length;
		double[][] y = neww;
		for (int row = 0; row < neww.length; row++) {
			for (int col = 0; col < neww[0].length; col++) {
				if (y[row][col] < ParRange_minnn[col]) {
					y[row][col] = 2 * ParRange_minnn[col] - y[row][col];
				}
				if (y[row][col] > ParRange_maxnn[col]) {
					y[row][col] = 2 * ParRange_maxnn[col] - y[row][col];
				}
			}
		}
		for (int row = 0; row < neww.length; row++) {
			for (int col = 0; col < neww[0].length; col++) {
				if (y[row][col] < ParRange_minn[col]) {
					y[row][col] = ParRange_minn[col] + Math.random()
							* (ParRange_maxn[col] - ParRange_minn[col]);
				}
				if (y[row][col] > ParRange_maxn[col]) {
					y[row][col] = ParRange_minn[col] + Math.random()
							* (ParRange_maxn[col] - ParRange_minn[col]);
				}
			}
		}

		return y;
	}

}