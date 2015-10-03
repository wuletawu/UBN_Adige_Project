
package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.jgrasstools.gears.libs.exceptions.ModelsIllegalargumentException;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.core.Dams;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.core.HillSlopeDuffy;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.core.Hydrometers;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.core.IDischargeContributor;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.core.IHillSlope;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.core.Offtakes;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.core.Tributaries;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.duffy.DuffyAdigeEngine;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.duffy.DuffyInputs;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.hymod.HymodAdigeEngineModified;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.utils.AdigeUtilities2;
import org.jgrasstools.hortonmachine.modules.network.PfafstetterNumber;
import org.joda.time.DateTime;
import org.joda.time.format.DateTimeFormatter;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

// Referenced classes of package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige:
//            IAdigeEngine

public class AdigeModified_2 extends JGTModel {

    public AdigeModified_2() {
        fNetnum = null;
        pNetNumCali = -1;
        fBaricenter = null;
        pRainintensity = -1D;
        pRainduration = -1;
        pPfafids = null;
        fMonpointid = null;
        fPfaff = null;
        fNetelevstart = null;
        fNetelevend = null;
        pCmax = null;
        pB = null;
        pAlpha = null;
        pRs = null;
        pDoReadMeas = false;
        pPathtoMeas = null;
        pMrain = Double.valueOf(1.0D);
        pMetp = Double.valueOf(1.0D);
        pRq = null;
        pQ0 = null;
        pLambda1 = null;
        pLambda2 = null;
        pQr = null;
        pBetanolinear = null;
        vR = null;
        pDoGeom = Boolean.valueOf(false);
        inSubbasinDist = null;
        fIdDistanze = null;
        fIdnetnum = null;
        pDorouting = false;
        pDimMeas = 0;
        doLog = false;
        tTimestep = 0;
        tStart = null;
        tEnd = null;
        inDuffyInput = null;
        pm = new LogProgressMonitor();
        rainArray = null;
        initialConditions = null;
        adigeEngine = null;
        outletHillslopeId = -1;
        conta = 0;
        init = true;
        init2 = true;
    }

    public void process()
        throws Exception {
        if(init) {
            outSimulated = new double[pDimMeas];
            outMeasured = new double[pDimMeas];
            init = false;
        }
        if(init2) {
            init2 = false;
            if(pDoReadMeas) {
                int dim = pDimMeas;
                double portate[] = new double[dim];
                int cont_portate = 0;
                try {
                    String str = new String();
                    str = pPathtoMeas;
                    FileInputStream fstream = new FileInputStream(str);
                    DataInputStream in = new DataInputStream(fstream);
                    BufferedReader br = new BufferedReader(new InputStreamReader(in));
                    double aa = 0.0D;
                    String strLine;
                    while((strLine = br.readLine()) != null)  {
                        aa = Double.parseDouble(strLine);
                        portate[cont_portate] = aa;
                        cont_portate++;
                    }
                    in.close();
                }
                catch(Exception exception) { }
                outMeasured = portate;
                pDoReadMeas = false;
            }
        }
        if(startTimestamp == null) {
            outDischarge = new HashMap();
            outSubdischarge = new HashMap();
            outVelocity = new HashMap();
            startTimestamp = adigeFormatter.parseDateTime(tStart);
            endTimestamp = adigeFormatter.parseDateTime(tEnd);
            currentTimstamp = startTimestamp;
            if(pRainintensity != -1D)
                if(pRainduration != -1)
                    rainEndTimestamp = startTimestamp.plusMinutes(pRainduration);
                else
                    throw new ModelsIllegalargumentException("In the case of usage of a constant rainintensity it is necessary to define also its duration.\nCheck your arguments, probably the --rainduration flag is missing.", this);
            if(fNetnum == null || fNetnum.length() < 1)
                throw new ModelsIllegalargumentException("Missing net num attribute name.", getClass().getSimpleName());
            if(fPfaff == null || fPfaff.length() < 1)
                throw new ModelsIllegalargumentException("Missing pfafstetter attribute name.", this);
            if(fMonpointid == null || fMonpointid.length() < 1)
                throw new ModelsIllegalargumentException("Missing monitoring point id attribute name.", getClass().getSimpleName());
            if(fBaricenter == null || fBaricenter.length() < 1)
                throw new ModelsIllegalargumentException("Missing basin centroid attribute name.", this);
            if(fNetelevstart == null || fNetelevstart.length() < 1)
                throw new ModelsIllegalargumentException("Missing start net elevation attribute name.", getClass().getSimpleName());
            if(fNetelevend == null || fNetelevend.length() < 1)
                throw new ModelsIllegalargumentException("Missing start net elevation attribute name.", getClass().getSimpleName());
            if(inHydrometers != null && inHydrometerdata != null && hydrometersHandler == null) {
                pm.message("Reading hydrometers geometries and mapping them to the network...");
                hydrometer_pfaff2idMap = new HashMap();
                hydrometersHandler = new Hydrometers(hydrometer_pfaff2idMap);
                FeatureIterator hydrometersIterator = inHydrometers.features();
                int pfaffIndex = -1;
                int monIdIndex = -1;
                String pNumberStr;
                int id;
                for(; hydrometersIterator.hasNext(); hydrometer_pfaff2idMap.put(pNumberStr, Integer.valueOf(id))) {
                    SimpleFeature hydrometer = (SimpleFeature)hydrometersIterator.next();
                    if(pfaffIndex == -1) {
                        SimpleFeatureType featureType = hydrometer.getFeatureType();
                        pfaffIndex = featureType.indexOf(fPfaff);
                        if(pfaffIndex == -1)
                            throw new ModelsIllegalargumentException((new StringBuilder("The hydrometer features are missing the pafaffstetter attribute field: ")).append(fPfaff).toString(), this);
                        monIdIndex = featureType.indexOf(fMonpointid);
                        if(monIdIndex == -1)
                            throw new ModelsIllegalargumentException((new StringBuilder("The hydrometer features are missing the id attribute field: ")).append(fMonpointid).toString(), this);
                    }
                    pNumberStr = hydrometer.getAttribute(pfaffIndex).toString();
                    id = ((Number)hydrometer.getAttribute(monIdIndex)).intValue();
                }

            }
            if(inDams != null && inDamsdata != null && damsHandler == null) {
                pm.message("Reading dams geometries and mapping them to the network...");
                dams_pfaff2idMap = new HashMap();
                damsHandler = new Dams(dams_pfaff2idMap);
                FeatureIterator damsIterator = inDams.features();
                int pfaffIndex = -1;
                int monIdIndex = -1;
                String pNumberStr;
                int id;
                for(; damsIterator.hasNext(); dams_pfaff2idMap.put(pNumberStr, Integer.valueOf(id))) {
                    SimpleFeature dam = (SimpleFeature)damsIterator.next();
                    if(pfaffIndex == -1) {
                        SimpleFeatureType featureType = dam.getFeatureType();
                        pfaffIndex = featureType.indexOf(fPfaff);
                        if(pfaffIndex == -1)
                            throw new ModelsIllegalargumentException((new StringBuilder("The dams features are missing the pfaffstetter attribute field: ")).append(fPfaff).toString(), this);
                        monIdIndex = featureType.indexOf(fMonpointid);
                        if(monIdIndex == -1)
                            throw new ModelsIllegalargumentException((new StringBuilder("The dams features are missing the id attribute field: ")).append(fMonpointid).toString(), this);
                    }
                    pNumberStr = (String)dam.getAttribute(pfaffIndex);
                    id = ((Number)dam.getAttribute(monIdIndex)).intValue();
                }

            }
            if(inTributary != null && inTributarydata != null && tributaryHandler == null) {
                pm.message("Reading tributary geometries and mapping them to the network...");
                tributary_pfaff2idMap = new HashMap();
                tributaryHandler = new Tributaries(tributary_pfaff2idMap);
                FeatureIterator tributaryIterator = inTributary.features();
                int pfaffIndex = -1;
                int monIdIndex = -1;
                String pNumberStr;
                int id;
                for(; tributaryIterator.hasNext(); tributary_pfaff2idMap.put(pNumberStr, Integer.valueOf(id))) {
                    SimpleFeature tributary = (SimpleFeature)tributaryIterator.next();
                    if(pfaffIndex == -1) {
                        SimpleFeatureType featureType = tributary.getFeatureType();
                        pfaffIndex = featureType.indexOf(fPfaff);
                        if(pfaffIndex == -1)
                            throw new ModelsIllegalargumentException((new StringBuilder("The tributary features are missing the pfaffstetter attribute field: ")).append(fPfaff).toString(), this);
                        monIdIndex = featureType.indexOf(fMonpointid);
                        if(monIdIndex == -1)
                            throw new ModelsIllegalargumentException((new StringBuilder("The tributary features are missing the id attribute field: ")).append(fMonpointid).toString(), this);
                    }
                    pNumberStr = (String)tributary.getAttribute(pfaffIndex);
                    id = ((Number)tributary.getAttribute(monIdIndex)).intValue();
                }

            }
            if(inOfftakes != null && inOfftakesdata != null && offtakesHandler == null) {
                pm.message("Reading offtakes geometries and mapping them to the network...");
                offtakes_pfaff2idMap = new HashMap();
                offtakesHandler = new Offtakes(offtakes_pfaff2idMap, pm);
                FeatureIterator offtakesIterator = inOfftakes.features();
                int pfaffIndex = -1;
                int monIdIndex = -1;
                String pNumberStr;
                int id;
                for(; offtakesIterator.hasNext(); offtakes_pfaff2idMap.put(pNumberStr, Integer.valueOf(id))) {
                    SimpleFeature offtakes = (SimpleFeature)offtakesIterator.next();
                    if(pfaffIndex == -1) {
                        SimpleFeatureType featureType = offtakes.getFeatureType();
                        pfaffIndex = featureType.indexOf(fPfaff);
                        if(pfaffIndex == -1)
                            throw new ModelsIllegalargumentException((new StringBuilder("The offtakes features are missing the pfaffstetter attribute field: ")).append(fPfaff).toString(), this);
                        monIdIndex = featureType.indexOf(fMonpointid);
                        if(monIdIndex == -1)
                            throw new ModelsIllegalargumentException((new StringBuilder("The offtakes features are missing the id attribute field: ")).append(fMonpointid).toString(), this);
                    }
                    pNumberStr = (String)offtakes.getAttribute(pfaffIndex);
                    id = ((Number)offtakes.getAttribute(monIdIndex)).intValue();
                }

            }
            hillsSlopeNum = inHillslope.size();
            orderedHillslopes = AdigeUtilities2.generateHillSlopes(inNetwork, inHillslope, fNetnum, fPfaff, fNetelevstart, fNetelevend, fBaricenter, pm);
            if(inDuffyInput != null) {
                List duffyHillslopes = new ArrayList();
                IHillSlope newHS;
                for(Iterator iterator = orderedHillslopes.iterator(); iterator.hasNext(); duffyHillslopes.add(newHS)) {
                    IHillSlope hillSlope = (IHillSlope)iterator.next();
                    newHS = new HillSlopeDuffy(hillSlope, inDuffyInput);
                }

                orderedHillslopes = duffyHillslopes;
            }
            IHillSlope outletHillSlope = (IHillSlope)orderedHillslopes.get(0);
            outletHillslopeId = outletHillSlope.getHillslopeId();
            netPfaffsList = new ArrayList();
            pfaff2Index = new HashMap();
            basinid2Index = new HashMap();
            index2Basinid = new HashMap();
            pm.beginTask("Analaysing hillslopes and calculating distribution curves...", orderedHillslopes.size());
            for(int i = 0; i < orderedHillslopes.size(); i++) {
                IHillSlope hillSlope = (IHillSlope)orderedHillslopes.get(i);
                PfafstetterNumber pfafstetterNumber = hillSlope.getPfafstetterNumber();
                netPfaffsList.add(pfafstetterNumber);
                int hillslopeId = hillSlope.getHillslopeId();
                basinid2Index.put(Integer.valueOf(hillslopeId), Integer.valueOf(i));
                index2Basinid.put(Integer.valueOf(i), Integer.valueOf(hillslopeId));
                pfaff2Index.put(pfafstetterNumber.toString(), Integer.valueOf(i));
                pm.worked(1);
            }

            pm.done();
            if(pPfafids == null)
                pPfafids = outletHillSlope.getPfafstetterNumber().toString();
            if(pfaffsList == null) {
                String split[] = pPfafids.split(",");
                for(int i = 0; i < split.length; i++)
                    split[i] = split[i].trim();

                pfaffsList = Arrays.asList(split);
            }
            if(pmode == 0) {
                initialConditions = new double[hillsSlopeNum * 4];
                adigeEngine = new DuffyAdigeEngine(orderedHillslopes, inDuffyInput, pm, doLog, initialConditions, basinid2Index, index2Basinid, pfaffsList, pfaff2Index, outDischarge, outSubdischarge, startTimestamp, endTimestamp, tTimestep);
            } else
            if(pmode == 1) {
                initialConditions = null;
                adigeEngine = new HymodAdigeEngineModified(pAlpha, pB, pCmax, pDoGeom.booleanValue(), pDorouting, pLambda1, pLambda2, pQr, fIdDistanze, vR, pRq, pRs, pMetp, pMrain, fIdnetnum, inSubbasinDist, pQ0, orderedHillslopes, index2Basinid, outDischarge, outSubdischarge, outVelocity, pfaffsList, doLog, doLog, pm, initialConditions);
            } else {
                throw new ModelsIllegalargumentException("No parameters for any model were defined. Check your syntax.", this);
            }
            if(hydrometersHandler != null)
                adigeEngine.addDischargeContributor(hydrometersHandler);
            if(damsHandler != null)
                adigeEngine.addDischargeContributor(damsHandler);
            if(tributaryHandler != null)
                adigeEngine.addDischargeContributor(tributaryHandler);
            if(offtakesHandler != null)
                adigeEngine.addDischargeContributor(offtakesHandler);
        } else {
            currentTimstamp = currentTimstamp.plusMinutes(tTimestep);
        }
        if(inHydrometerdata != null)
            hydrometersHandler.setCurrentData(inHydrometerdata);
        if(inDamsdata != null)
            damsHandler.setCurrentData(inDamsdata);
        if(inOfftakesdata != null)
            offtakesHandler.setCurrentData(inOfftakesdata);
        if(inTributarydata != null)
            tributaryHandler.setCurrentData(inTributarydata);
        if(pRainintensity != -1D) {
            rainArray = new double[netPfaffsList.size()];
            if(currentTimstamp.isBefore(rainEndTimestamp))
                Arrays.fill(rainArray, pRainintensity);
            else
                Arrays.fill(rainArray, 0.0D);
        } else {
            rainArray = new double[hillsSlopeNum];
            etpArray = new double[hillsSlopeNum];
            setDataArray(inRain, rainArray);
            if(inEtp != null)
                setDataArray(inEtp, etpArray);
        }
        
    //     long runningDateInMinutes = currentTimstamp.getMillis() / 1000L / 60L;
    //     double intervalStartTimeInMinutes = runningDateInMinutes;
    //     double intervalEndTimeInMinutes = runningDateInMinutes + tTimestep;

        initialConditions = adigeEngine.solve(currentTimstamp, tTimestep, 1.0D, initialConditions, rainArray, etpArray);
        if(pNetNumCali != -1 && conta < outSimulated.length) {
            outSimulated[conta] = ((double[])outDischarge.get(Integer.valueOf(pNetNumCali)))[0];
            conta++;
        }
    }

    private void setDataArray(HashMap dataMap, double endArray[]) {
        Set entries = dataMap.entrySet();
        for(Iterator iterator = entries.iterator(); iterator.hasNext();) {
            java.util.Map.Entry entry = (java.util.Map.Entry)iterator.next();
            Integer id = (Integer)entry.getKey();
            double value[] = (double[])entry.getValue();
            Integer index = (Integer)basinid2Index.get(id);
            if(index != null) {
                if(JGTConstants.isNovalue(value[0]))
                    value[0] = 0.0D;
                endArray[index.intValue()] = value[0];
            }
        }

    }

    public SimpleFeatureCollection inHillslope;
    public String fNetnum;
    public int pNetNumCali;
    public String fBaricenter;
    public double pRainintensity;
    public int pRainduration;
    public HashMap inRain;
    public SimpleFeatureCollection inHydrometers;
    public HashMap inHydrometerdata;
    public SimpleFeatureCollection inDams;
    public HashMap inDamsdata;
    public SimpleFeatureCollection inTributary;
    public HashMap inTributarydata;
    public SimpleFeatureCollection inOfftakes;
    public HashMap inOfftakesdata;
    public String pPfafids;
    public String fMonpointid;
    public SimpleFeatureCollection inNetwork;
    public String fPfaff;
    public int pmode;
    public String fNetelevstart;
    public String fNetelevend;
    public Double pCmax;
    public Double pB;
    public Double pAlpha;
    public Double pRs;
    public boolean pDoReadMeas;
    public String pPathtoMeas;
    public Double pMrain;
    public Double pMetp;
    public Double pRq;
    public Double pQ0;
    public Double pLambda1;
    public Double pLambda2;
    public Double pQr;
    public Double pBetanolinear;
    public Double vR;
    public Boolean pDoGeom;
    public SimpleFeatureCollection inSubbasinDist;
    public String fIdDistanze;
    public String fIdnetnum;
    public boolean pDorouting;
    public int pDimMeas;
    public HashMap inEtp;
    public boolean doLog;
    public int tTimestep;
    public String tStart;
    public String tEnd;
    public DuffyInputs inDuffyInput;
    public IJGTProgressMonitor pm;
    public HashMap outDischarge;
    public HashMap outSubdischarge;
    public HashMap outVelocity;
    public double outMeasured[];
    public double outSimulated[];
    private double rainArray[];
    private double etpArray[];
    private double initialConditions[];
    private IAdigeEngine adigeEngine;
    private List netPfaffsList;
    private IDischargeContributor hydrometersHandler;
    private HashMap hydrometer_pfaff2idMap;
    private IDischargeContributor damsHandler;
    private HashMap dams_pfaff2idMap;
    private IDischargeContributor tributaryHandler;
    private HashMap tributary_pfaff2idMap;
    private IDischargeContributor offtakesHandler;
    private HashMap offtakes_pfaff2idMap;
    private HashMap basinid2Index;
    private HashMap index2Basinid;
    private int hillsSlopeNum;
    private int outletHillslopeId;
    private HashMap pfaff2Index;
    private List orderedHillslopes;
    private int conta;
    public static DateTimeFormatter adigeFormatter;
    private DateTime startTimestamp;
    private DateTime endTimestamp;
    private DateTime currentTimstamp;
    private DateTime rainEndTimestamp;
    private List pfaffsList;
    boolean init;
    boolean init2;

    static  {
        adigeFormatter = JGTConstants.utcDateFormatterYYYYMMDDHHMM;
    }
}
