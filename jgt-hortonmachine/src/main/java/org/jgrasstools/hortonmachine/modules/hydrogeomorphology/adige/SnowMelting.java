/*package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige;

import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Set;

import javax.media.jai.iterator.RandomIterFactory;
import javax.media.jai.iterator.WritableRandomIter;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.geometry.DirectPosition2D;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;
import org.jgrasstools.gears.utils.CrsUtilities;
import org.jgrasstools.gears.utils.RegionMap;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.gears.utils.geometry.GeometryUtilities;
import org.joda.time.DateTime;
import org.joda.time.DateTimeZone;
import org.joda.time.format.DateTimeFormat;
import org.joda.time.format.DateTimeFormatter;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.geometry.DirectPosition;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Point;

public class SnowMelting extends JGTModel {

    public SnowMelting() {
        inSkyview = null;
        pm1 = 1.0D;
        tStartDate = null;
        tEndDate = null;
        inDem = null;
        time = null;
        inTempGrid = null;
        incondIniIGrid = null;
        incondIniLGrid = null;
        inRainGrid = null;
        inInsFeb = null;
        inInsJan = null;
        inInsMar = null;
        inInsApr = null;
        inInsMay = null;
        inInsJun = null;
        inStations = null;
        fStationsid = null;
        doRaster = true;
        numSitesWhereCalibrate = 1;
        pDoReadMeas = false;
        pm = new LogProgressMonitor();
        outGrid = null;
        pathToSwe = null;
        pathToPrainOutput = null;
        pathToSolidWater = null;
        pathToLiquidWater = null;
        pPathtoMeas = null;
        outMeasured = null;
        pathToMelting = null;
        outMeltingDataGrid = null;
        outSweDataGrid = null;
        outSwevector = null;
        height = 0;
        width = 0;
        outMeltingData = null;
        outSWEData = null;
        outPrainData = null;
        insFebWR = null;
        insJanWR = null;
        intempWR = null;
        inrainWR = null;
        demwr = null;
        insMarWR = null;
        insAprWR = null;
        insJunWR = null;
        insMayWR = null;
        flag = true;
        skyviewfactorWR = null;
        init2 = true;
        hoursToJan = 720D;
        hoursToFeb = 1416D;
        hoursToMar = 2160D;
        hoursToApr = 2880D;
        hoursToMay = 3631D;
        hoursToJune = 4351D;
        outWRSWE = null;
        outWRMEL = null;
        outMELTINGVECTOR = null;
        outSWEVECTOR = null;
        outMELTINGVECTOR_forMaps = null;
        outSWEVECTOR_forMaps = null;
        folder = null;
        P_listOfFiles = null;
        conta = 0;
        pointsToComputeId2Coordinates = null;
        doOneTime = true;
    }

    public void process()
        throws Exception {
        if(init2) {
            outMeasured = new double[pDimMeas];
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
        CoordinateReferenceSystem sourceCRS = inDem.getCoordinateReferenceSystem2D();
        int numPointToCompute = 0;
        GridGeometry2D inRainGridGeo = null;
        if(!doRaster) {
            pointsToComputeId2Coordinates = getCoordinate(numPointToCompute, inStations, fStationsid);
            numPointToCompute = pointsToComputeId2Coordinates.size();
        } else
        if(doRaster) {
            if(inRainGrid != null)
                inRainGridGeo = inRainGrid.getGridGeometry();
            pointsToComputeId2Coordinates = getCoordinate(inRainGridGeo);
            numPointToCompute = pointsToComputeId2Coordinates.size();
        }
        Set pointsToInterpolateIdSet = pointsToComputeId2Coordinates.keySet();
        Iterator idIterator = pointsToInterpolateIdSet.iterator();
        int j = 0;
        xStation = new double[numPointToCompute];
        yStation = new double[numPointToCompute];
        idStation = new int[numPointToCompute];
        colnumvetVect = new int[numPointToCompute];
        rownumvetVect = new int[numPointToCompute];
        lambdaVect = new double[numPointToCompute];
        CoordinateReferenceSystem targetCRS = DefaultGeographicCRS.WGS84;
        while(idIterator.hasNext())  {
            int id = ((Integer)idIterator.next()).intValue();
            idStation[j] = id;
            Coordinate coordinate = (Coordinate)pointsToComputeId2Coordinates.get(Integer.valueOf(id));
            xStation[j] = coordinate.x;
            yStation[j] = coordinate.y;
            double srcPts[] = {
                xStation[j], yStation[j]
            };
            Coordinate source = new Coordinate(srcPts[0], srcPts[1]);
            Point so[] = {
                GeometryUtilities.gf().createPoint(source)
            };
            CrsUtilities.reproject(sourceCRS, targetCRS, so);
            lambdaVect[j] = Math.toRadians(so[0].getY());
            j++;
        }
        MathTransform transf = inDem.getGridGeometry().getCRSToGrid2D();
        for(int i = 0; i < xStation.length; i++) {
            DirectPosition point = new DirectPosition2D(sourceCRS, xStation[i], yStation[i]);
            DirectPosition gridPoint = transf.transform(point, null);
            colnumvetVect[i] = (int)gridPoint.getCoordinate()[0];
            rownumvetVect[i] = (int)gridPoint.getCoordinate()[1];
        }

        double minimofeb = 0.0D;
        double minimomar = 0.0D;
        double minimoapr = 0.0D;
        double minimomagg = 0.0D;
        double minimogiu = 0.0D;
        double dx = 0.0D;
        DateTimeFormatter formatter = DateTimeFormat.forPattern("yyyy-MM-dd HH:mm").withZone(DateTimeZone.UTC);
        DateTime startcurrentDatetime = formatter.parseDateTime(tStartDate);
        DateTime endcurrentDatetime = formatter.parseDateTime(tEndDate);
        long diff = 0L;
        if(!doDaily)
            diff = (endcurrentDatetime.getMillis() - startcurrentDatetime.getMillis()) / 0x36ee80L;
        if(doDaily)
            diff = (endcurrentDatetime.getMillis() - startcurrentDatetime.getMillis()) / 0x5265c00L;
        DateTime array[] = new DateTime[(int)diff];
        if(!doDaily) {
            for(int i = 0; i < array.length; i++)
                array[i] = startcurrentDatetime.plusHours(i);

        }
        if(doDaily) {
            for(int i = 0; i < array.length; i++)
                array[i] = startcurrentDatetime.plusDays(i);

        }
        if(doOneTime) {
            attribute = CoverageUtilities.getRegionParamsFromGridCoverage(inSkyview);
            dx = ((Double)attribute.get("XRES")).doubleValue();
            double srcPts[] = {
                ((Double)attribute.get("EAST")).doubleValue(), ((Double)attribute.get("SOUTH")).doubleValue()
            };
            Coordinate source = new Coordinate(srcPts[0], srcPts[1]);
            Point so[] = {
                GeometryUtilities.gf().createPoint(source)
            };
            CrsUtilities.reproject(sourceCRS, targetCRS, so);
            lambda = Math.toRadians(so[0].getY());
            CoverageUtilities.getRegionParamsFromGridCoverage(inSkyview);
            if(pMode == 0) {
                RenderedImage insFebTmpRI = inInsFeb.getRenderedImage();
                width = insFebTmpRI.getWidth();
                height = insFebTmpRI.getHeight();
                insFebWR = CoverageUtilities.replaceNovalue(insFebTmpRI, -9999D);
                RenderedImage insJanTmpRI = inInsJan.getRenderedImage();
                width = insJanTmpRI.getWidth();
                height = insJanTmpRI.getHeight();
                insJanWR = CoverageUtilities.replaceNovalue(insJanTmpRI, -9999D);
                RenderedImage insMarTmpRI = inInsMar.getRenderedImage();
                insMarWR = CoverageUtilities.replaceNovalue(insMarTmpRI, -9999D);
                insMarTmpRI = null;
                RenderedImage insAprTmpRI = inInsApr.getRenderedImage();
                insAprWR = CoverageUtilities.replaceNovalue(insAprTmpRI, -9999D);
                insAprTmpRI = null;
                RenderedImage insMayTmpRI = inInsMay.getRenderedImage();
                insMayWR = CoverageUtilities.replaceNovalue(insMayTmpRI, -9999D);
                insMayTmpRI = null;
                RenderedImage insJuneTmpRI = inInsJun.getRenderedImage();
                insJunWR = CoverageUtilities.replaceNovalue(insJuneTmpRI, -9999D);
                insJuneTmpRI = null;
                minimofeb = findminumum(insFebWR, hoursToFeb);
                minimomar = findminumum(insMarWR, hoursToMar);
                minimoapr = findminumum(insAprWR, hoursToApr);
                minimomagg = findminumum(insMayWR, hoursToMay);
                minimogiu = findminumum(insJunWR, hoursToJune);
            }
            RenderedImage drmri = inDem.getRenderedImage();
            demwr = CoverageUtilities.replaceNovalue(drmri, -9999D);
            drmri = null;
            RenderedImage SkyTmpRI = inSkyview.getRenderedImage();
            skyviewfactorWR = CoverageUtilities.renderedImage2WritableRaster(SkyTmpRI, true);
            inCondIniL = new HashMap();
            inCondIniI = new HashMap();
            outSWEData = new HashMap();
            outPrainData = new HashMap();
            doOneTime = false;
        }
        folder = null;
        P_listOfFiles = null;
        T_listOfFiles = null;
        if(doRaster) {
            folder = new File(pathRainfMaps);
            P_listOfFiles = folder.listFiles();
            for(int i = 0; i < P_listOfFiles.length; i++)
                System.out.println(P_listOfFiles[i]);

            folder = new File(pathTempMaps);
            T_listOfFiles = folder.listFiles();
            for(int i = 0; i < T_listOfFiles.length; i++)
                System.out.println(T_listOfFiles[i]);

        }
        TimeSeriesIteratorReader reader_temp = new TimeSeriesIteratorReader();
        TimeSeriesIteratorReader reader_solar = new TimeSeriesIteratorReader();
        TimeSeriesIteratorReader reader_rainf = new TimeSeriesIteratorReader();
        TimeSeriesIteratorWriter writer = new TimeSeriesIteratorWriter();
        TimeSeriesIteratorWriter writer2 = new TimeSeriesIteratorWriter();
        TimeSeriesIteratorWriter writer3 = new TimeSeriesIteratorWriter();
        TimeSeriesIteratorWriter writer4 = new TimeSeriesIteratorWriter();
        TimeSeriesIteratorWriter writer5 = new TimeSeriesIteratorWriter();
        if(!doRaster) {
            if(pathTemp != null) {
                reader_temp.file = pathTemp;
                reader_temp.idfield = "ID";
                reader_temp.tStart = tStartDate;
                reader_temp.tEnd = tEndDate;
                reader_temp.fileNovalue = "-9999";
                reader_temp.tTimestep = inTimestep;
            }
            if(pathToSolarRad != null) {
                reader_solar.file = pathToSolarRad;
                reader_solar.idfield = "ID";
                reader_solar.tStart = tStartDate;
                reader_solar.tEnd = tEndDate;
                reader_solar.fileNovalue = "-9999";
                reader_solar.tTimestep = inTimestep;
            }
            if(pathRainf != null) {
                reader_rainf.file = pathRainf;
                reader_rainf.idfield = "ID";
                reader_rainf.tStart = tStartDate;
                reader_rainf.tEnd = tEndDate;
                reader_rainf.fileNovalue = "-9999";
                reader_rainf.tTimestep = inTimestep;
            }
            if(pathToMelting != null) {
                writer.file = pathToMelting;
                writer.tStart = tStartDate;
                writer.tTimestep = inTimestep;
            }
            if(pathToSwe != null) {
                writer2.file = pathToSwe;
                writer2.tStart = tStartDate;
                writer2.tTimestep = inTimestep;
            }
            if(pathToPrainOutput != null) {
                writer3.file = pathToPrainOutput;
                writer3.tStart = tStartDate;
                writer3.tTimestep = inTimestep;
            }
            if(pathToSolidWater != null) {
                writer4.file = pathToSolidWater;
                writer4.tStart = tStartDate;
                writer4.tTimestep = inTimestep;
            }
            if(pathToLiquidWater != null) {
                writer5.file = pathToLiquidWater;
                writer5.tStart = tStartDate;
                writer5.tTimestep = inTimestep;
            }
            if(numSitesWhereCalibrate != 1) {
                outSWEVECTOR = new double[array.length * numSitesWhereCalibrate];
                outMELTINGVECTOR = new double[array.length * numSitesWhereCalibrate];
            } else {
                outSWEVECTOR = new double[array.length];
                outMELTINGVECTOR = new double[array.length];
            }
            conta = 0;
        }
        outMELTINGVECTOR_forMaps = new double[numPointToCompute];
        outSWEVECTOR_forMaps = new double[numPointToCompute];
        for(int i = 0; i < array.length; i++) {
            outSWE = new HashMap();
            outMeltingData = new HashMap();
            DateTime currentime = array[i];
            if(!doRaster) {
                if(pathTemp != null) {
                    reader_temp.nextRecord();
                    temp_values = reader_temp.outData;
                }
                if(pathRainf != null) {
                    reader_rainf.nextRecord();
                    rain_values = reader_rainf.outData;
                }
                if(pathToSolarRad != null) {
                    reader_solar.nextRecord();
                    solar_values = reader_solar.outData;
                }
            }
            int month = currentime.getMonthOfYear();
            switch(month) {
            case 1: // '\001'
                calcMelting(insJanWR, demwr, dx, currentime, hoursToJan, minimofeb, i);
                break;

            case 10: // '\n'
                calcMelting(insFebWR, demwr, dx, currentime, hoursToFeb, minimofeb, i);
                break;

            case 11: // '\013'
                calcMelting(insFebWR, demwr, dx, currentime, hoursToFeb, minimofeb, i);
                break;

            case 12: // '\f'
                calcMelting(insFebWR, demwr, dx, currentime, hoursToFeb, minimofeb, i);
                break;

            case 2: // '\002'
                calcMelting(insFebWR, demwr, dx, currentime, hoursToFeb, minimofeb, i);
                break;

            case 3: // '\003'
                calcMelting(insMarWR, demwr, dx, currentime, hoursToMar, minimomar, i);
                break;

            case 4: // '\004'
                calcMelting(insAprWR, demwr, dx, currentime, hoursToApr, minimoapr, i);
                break;

            case 5: // '\005'
                calcMelting(insMayWR, demwr, dx, currentime, hoursToMay, minimomagg, i);
                break;

            case 6: // '\006'
                calcMelting(insJunWR, demwr, dx, currentime, hoursToJune, minimogiu, i);
                break;

            case 7: // '\007'
                calcMelting(insJunWR, demwr, dx, currentime, hoursToJune, minimogiu, i);
                break;

            case 8: // '\b'
                calcMelting(insJunWR, demwr, dx, currentime, hoursToJune, minimogiu, i);
                break;

            case 9: // '\t'
                calcMelting(insJunWR, demwr, dx, currentime, hoursToJune, minimogiu, i);
                break;
            }
            if(!doRaster) {
                if(pathToMelting != null) {
                    writer.inData = outMeltingData;
                    writer.writeNextLine();
                }
                if(pathToSolidWater != null) {
                    writer4.inData = inCondIniI;
                    writer4.writeNextLine();
                }
                if(pathToLiquidWater != null) {
                    writer5.inData = inCondIniL;
                    writer5.writeNextLine();
                }
                if(pathToSwe != null) {
                    writer2.inData = outSWEData;
                    writer2.writeNextLine();
                }
                if(pathToPrainOutput != null) {
                    writer3.inData = outPrainData;
                    writer3.writeNextLine();
                }
            } else
            if(doRaster)
                storeResult(outSWEVECTOR_forMaps, outMELTINGVECTOR_forMaps, pointsToComputeId2Coordinates);
        }

        if(pathToMelting != null)
            writer.close();
        if(pathToSwe != null)
            writer2.close();
        if(pathToPrainOutput != null)
            writer3.close();
        if(pathToSolidWater != null)
            writer4.close();
        if(pathToLiquidWater != null)
            writer5.close();
        outSwevector = outSWEVECTOR;
    }

    private void calcMelting(WritableRaster energyWR, WritableRaster demWR, double dx, DateTime time, double ore, 
            double minimoEI, int iii)
        throws Exception {
        CoverageUtilities.getRegionParamsFromGridCoverage(inSkyview);
        double condiniI = (0.0D / 0.0D);
        double condiniL = (0.0D / 0.0D);
        for(int contastaz = 0; contastaz < xStation.length; contastaz++) {
            int colnuber = colnumvetVect[contastaz];
            int rownumber = rownumvetVect[contastaz];
            int i = colnuber;
            int j = rownumber;
            int id = idStation[contastaz];
            double temperatura = (0.0D / 0.0D);
            double pioggia = (0.0D / 0.0D);
            double solar = (0.0D / 0.0D);
            if(!doRaster) {
                if(temp_values != null)
                    temperatura = ((double[])temp_values.get(Integer.valueOf(id)))[0];
                if(rain_values != null)
                    pioggia = ((double[])rain_values.get(Integer.valueOf(id)))[0];
                if(solar_values != null)
                    solar = ((double[])solar_values.get(Integer.valueOf(id)))[0];
                if(flag) {
                    if(condiniI_values == null)
                        condiniI = 0.0D;
                    if(condiniL_values == null)
                        condiniL = 0.0D;
                    if(contastaz == xStation.length - 1)
                        flag = false;
                } else {
                    condiniI = ((double[])inCondIniI.get(Integer.valueOf(id)))[0];
                    condiniL = ((double[])inCondIniL.get(Integer.valueOf(id)))[0];
                }
            } else
            if(doRaster) {
                RasterReader readersT = new RasterReader();
                readersT.file = T_listOfFiles[iii].toString();
                readersT.fileNovalue = Double.valueOf(-9999D);
                readersT.geodataNovalue = Double.valueOf((0.0D / 0.0D));
                readersT.process();
                GridCoverage2D t = readersT.outRaster;
                RasterReader readersP = new RasterReader();
                readersP.file = P_listOfFiles[iii].toString();
                readersP.fileNovalue = Double.valueOf(-9999D);
                readersP.geodataNovalue = Double.valueOf((0.0D / 0.0D));
                readersP.process();
                GridCoverage2D p = readersP.outRaster;
                RenderedImage inRainRI = p.getRenderedImage();
                inrainWR = CoverageUtilities.replaceNovalue(inRainRI, -9999D);
                inRainRI = null;
                RenderedImage inTempRI = t.getRenderedImage();
                intempWR = CoverageUtilities.replaceNovalue(inTempRI, -9999D);
                inTempRI = null;
                if(temp_values == null)
                    temperatura = intempWR.getSampleDouble(i, j, 0);
                if(rain_values == null)
                    pioggia = inrainWR.getSampleDouble(i, j, 0);
                if(flag) {
                    if(condiniI_values == null)
                        condiniI = 0.0D;
                    if(condiniL_values == null)
                        condiniL = 0.0D;
                    if(contastaz == xStation.length - 1)
                        flag = false;
                } else {
                    condiniI = ((double[])inCondIniI.get(Integer.valueOf(id)))[0];
                    condiniL = ((double[])inCondIniL.get(Integer.valueOf(id)))[0];
                }
            }
            double aaa[] = calcMelting(i, j, energyWR, demWR, temperatura, pioggia, solar, ore, time, minimoEI, condiniI, condiniL, lambdaVect[contastaz]);
            outSWEVECTOR_forMaps[contastaz] = aaa[0];
            inCondIniI.put(Integer.valueOf(id), new double[] {
                aaa[1]
            });
            inCondIniL.put(Integer.valueOf(id), new double[] {
                aaa[2]
            });
            outMeltingData.put(Integer.valueOf(id), new double[] {
                aaa[3]
            });
            outSWEData.put(Integer.valueOf(id), new double[] {
                aaa[0]
            });
            outPrainData.put(Integer.valueOf(id), new double[] {
                aaa[4]
            });
            outMELTINGVECTOR_forMaps[contastaz] = aaa[3];
            conta++;
        }

    }

    private double[] calcMelting(int i, int j, WritableRaster enWR, WritableRaster demWR, double tem, double pio, double sol, double ore, DateTime time, double minimoEI, double condiniI, double condiniL, double lll)
        throws IOException {
        double risultato[] = new double[5];
        double melting = 0.0D;
        double prain = 0.0D;
        double psnow = 0.0D;
        int day = time.getDayOfYear();
        double dayangb = 0.98562628336755642D * ((double)day - 79.436000000000007D);
        dayangb = Math.toRadians(dayangb);
        double delta = getDeclination(dayangb);
        double ss = Math.acos(-Math.tan(delta) * Math.tan(lll));
        double sunrise = 12D * (1.0D - ss / 3.1415926535897931D);
        double sunset = 12D * (1.0D + ss / 3.1415926535897931D);
        double hhh = (double)time.getMillisOfDay() / 3600000D;
        if(!doDaily) {
            if(hhh >= sunrise && hhh <= sunset) {
                if(skyviewfactorWR.getSampleDouble(i, j, 0) != -9999D) {
                    double ei = 0.0D;
                    if(pMode == 0)
                        ei = enWR.getSampleDouble(i, j, 0) / (ore / 24D);
                    double temp = tem;
                    if(JGTConstants.isNovalue(tem)) {
                        double zzz = demwr.getSampleDouble(i, j, 0);
                        temp = 273D + -0.0064999999999999997D * (zzz - 4000D);
                    }
                    if(JGTConstants.isNovalue(pio) || pio < 0.0D)
                        pio = 0.0D;
                    if(pMode == 0) {
                        if(temp > pTmelt) {
                            melting = ei * pCmf * (temp - pTmelt) * skyviewfactorWR.getSampleDouble(i, j, 0);
                            melting = Math.min(melting, condiniI + condiniL);
                        } else {
                            melting = 0.0D;
                        }
                    } else
                    if(pMode == 1) {
                        if(temp > pTmelt)
                            melting = pCmf * (temp - pTmelt) * skyviewfactorWR.getSampleDouble(i, j, 0);
                        else
                            melting = 0.0D;
                        melting = Math.max(melting, 0.0D);
                    } else
                    if(pMode == 2) {
                        if(temp > pTmelt)
                            melting = (pCmf + sol * pCrf) * (temp - pTmelt) * skyviewfactorWR.getSampleDouble(i, j, 0);
                        else
                            melting = 0.0D;
                        melting = Math.max(melting, 0.0D);
                    } else {
                        System.out.println("ERROREEEE");
                    }
                    prain = (pio / 3.1415926535897931D) * Math.atan((temp - pTmelt) / pm1) + pio * 0.5D;
                    psnow = pio - prain;
                    prain *= pCr;
                    psnow = pCs * psnow;
                    double freezing = 0.0D;
                    if(temp < pTmelt)
                        freezing = pCff * (pTmelt - temp);
                    else
                        freezing = 0.0D;
                    double I = 0.0D;
                    double deltat = 1.0D;
                    I = condiniI + deltat * ((psnow + freezing) - melting);
                    if(I < 0.0D) {
                        I = 0.0D;
                        melting = 0.0D;
                    }
                    double L = 0.0D;
                    L = condiniL + deltat * ((prain + melting) - freezing);
                    double L_max = pR * I;
                    double melting_discharge = 0.0D;
                    if(L < 0.0D)
                        L = 0.0D;
                    if(L > L_max) {
                        melting_discharge = L - L_max;
                        L = L_max;
                    }
                    risultato[0] = L + I;
                    risultato[1] = I;
                    risultato[2] = L;
                    risultato[3] = melting_discharge;
                    risultato[4] = freezing;
                    risultato[4] = prain;
                } else {
                    risultato[0] = -9999D;
                    risultato[1] = condiniI;
                    risultato[2] = condiniL;
                    risultato[3] = -9999D;
                    risultato[4] = -9999D;
                    risultato[4] = -9999D;
                }
            } else
            if(skyviewfactorWR.getSampleDouble(i, j, 0) != -9999D) {
                double z = 0.0D;
                if(pMode == 0)
                    z = minimoEI / ore;
                double temp = tem;
                if(JGTConstants.isNovalue(tem)) {
                    double zzz = demWR.getSampleDouble(i, j, 0);
                    temp = -0.0064999999999999997D * (zzz - 4000D);
                }
                if(pMode == 0) {
                    if(temp > pTmelt) {
                        melting = z * pCmf * (temp - pTmelt) * skyviewfactorWR.getSampleDouble(i, j, 0);
                        melting = Math.min(melting, condiniI + condiniL);
                    } else {
                        melting = 0.0D;
                    }
                } else
                if(pMode == 1) {
                    if(temp > pTmelt)
                        melting = pCmf * (temp - pTmelt) * skyviewfactorWR.getSampleDouble(i, j, 0);
                    else
                        melting = 0.0D;
                    melting = Math.max(melting, 0.0D);
                } else
                if(pMode == 2) {
                    if(temp > pTmelt)
                        melting = (pCmf + sol * pCrf) * (temp - pTmelt) * skyviewfactorWR.getSampleDouble(i, j, 0);
                    else
                        melting = 0.0D;
                    melting = Math.max(melting, 0.0D);
                } else {
                    System.out.println("ERROREEEE");
                }
                if(JGTConstants.isNovalue(pio) || pio < 0.0D)
                    pio = 0.0D;
                prain = (pio / 3.1415926535897931D) * Math.atan((temp - pTmelt) / pm1) + pio * 0.5D;
                psnow = pio - prain;
                prain *= pCr;
                psnow = pCs * psnow;
                double freezing = 0.0D;
                if(temp < pTmelt)
                    freezing = pCff * (pTmelt - temp);
                else
                    freezing = 0.0D;
                double I = 0.0D;
                double deltat = 1.0D;
                I = condiniI + deltat * ((psnow + freezing) - melting);
                if(I < 0.0D) {
                    I = 0.0D;
                    melting = 0.0D;
                }
                double L = 0.0D;
                L = condiniL + deltat * ((prain + melting) - freezing);
                double L_max = pR * I;
                double melting_discharge = 0.0D;
                if(L < 0.0D)
                    L = 0.0D;
                if(L > L_max) {
                    melting_discharge = L - L_max;
                    L = L_max;
                }
                risultato[0] = L + I;
                risultato[1] = I;
                risultato[2] = L;
                risultato[3] = melting_discharge;
                risultato[4] = freezing;
                risultato[4] = prain;
            } else {
                risultato[0] = -9999D;
                risultato[1] = condiniI;
                risultato[2] = condiniL;
                risultato[3] = -9999D;
                risultato[4] = -9999D;
                risultato[4] = -9999D;
            }
        } else
        if(skyviewfactorWR.getSampleDouble(i, j, 0) != -9999D) {
            double ei = 0.0D;
            if(pMode == 0)
                ei = enWR.getSampleDouble(i, j, 0) / (ore / 24D);
            double temp = tem;
            if(JGTConstants.isNovalue(tem)) {
                double zzz = demWR.getSampleDouble(i, j, 0);
                temp = 273D + -0.0064999999999999997D * (zzz - 4000D);
            }
            if(JGTConstants.isNovalue(pio))
                pio = 0.0D;
            if(pMode == 0) {
                if(temp > pTmelt)
                    melting = ei * pCmf * (temp - pTmelt) * skyviewfactorWR.getSampleDouble(i, j, 0);
                else
                    melting = 0.0D;
                melting = Math.max(melting, 0.0D);
            } else
            if(pMode == 1) {
                if(temp > pTmelt)
                    melting = pCmf * (temp - pTmelt) * skyviewfactorWR.getSampleDouble(i, j, 0);
                else
                    melting = 0.0D;
                melting = Math.max(melting, 0.0D);
            } else
            if(pMode == 2) {
                if(temp > pTmelt)
                    melting = (pCmf + sol * pCrf) * (temp - pTmelt) * skyviewfactorWR.getSampleDouble(i, j, 0);
                else
                    melting = 0.0D;
                melting = Math.max(melting, 0.0D);
            } else {
                System.out.println("ERROREEEE");
            }
            prain = (pio / 3.1415926535897931D) * Math.atan((temp - pTmelt) / pm1) + pio * 0.5D;
            psnow = pio - prain;
            prain *= pCr;
            psnow = pCs * psnow;
            double freezing = 0.0D;
            if(temp < pTmelt)
                freezing = pCff * (pTmelt - temp);
            else
                freezing = 0.0D;
            double I = 0.0D;
            double deltat = 1.0D;
            I = condiniI + deltat * ((psnow + freezing) - melting);
            if(I < 0.0D) {
                I = 0.0D;
                melting = 0.0D;
            }
            double L = 0.0D;
            L = condiniL + deltat * ((prain + melting) - freezing);
            double L_max = pR * I;
            double melting_discharge = 0.0D;
            if(L < 0.0D)
                L = 0.0D;
            if(L > L_max) {
                melting_discharge = L - L_max;
                L = L_max;
            }
            risultato[0] = L + I;
            risultato[1] = I;
            risultato[2] = L;
            risultato[3] = melting_discharge;
            risultato[4] = freezing;
            risultato[4] = prain;
        } else {
            risultato[0] = -9999D;
            risultato[1] = condiniI;
            risultato[2] = condiniL;
            risultato[3] = -9999D;
            risultato[4] = -9999D;
            risultato[4] = -9999D;
        }
        return risultato;
    }

    private double findminumum(WritableRaster w, double num) {
        double minimo = 10000000D;
        for(int j = 0; j < height; j++) {
            for(int i = 0; i < width; i++)
                if(w.getSample(i, j, 0) != -9999 && (double)w.getSample(i, j, 0) / num > 0.0D && (double)w.getSample(i, j, 0) / num < minimo)
                    minimo = (double)w.getSample(i, j, 0) / num;

        }

        return minimo;
    }

    private LinkedHashMap getCoordinate(GridGeometry2D grid) {
        LinkedHashMap out = new LinkedHashMap();
        int count = 0;
        RegionMap regionMap = CoverageUtilities.gridGeometry2RegionParamsMap(grid);
        cols = regionMap.getCols();
        rows = regionMap.getRows();
        south = regionMap.getSouth();
        west = regionMap.getWest();
        xres = regionMap.getXres();
        yres = regionMap.getYres();
        outWRSWE = CoverageUtilities.createDoubleWritableRaster(cols, rows, null, null, null);
        outWRMEL = CoverageUtilities.createDoubleWritableRaster(cols, rows, null, null, null);
        double northing = south;
        double easting = west;
        for(int i = 0; i < cols; i++) {
            easting += xres;
            for(int j = 0; j < rows; j++) {
                northing += yres;
                Coordinate coordinate = new Coordinate();
                coordinate.x = west + (double)i * xres;
                coordinate.y = south + (double)j * yres;
                out.put(Integer.valueOf(count), coordinate);
                count++;
            }

        }

        return out;
    }

    private LinkedHashMap getCoordinate(int nStaz, SimpleFeatureCollection collection, String idField)
        throws Exception {
        LinkedHashMap id2CoordinatesMap;
        FeatureIterator iterator;
        id2CoordinatesMap = new LinkedHashMap();
        iterator = collection.features();
        Coordinate coordinate = null;
        Coordinate coordinate;
        int name;
        for(; iterator.hasNext(); id2CoordinatesMap.put(Integer.valueOf(name), coordinate)) {
            SimpleFeature feature = (SimpleFeature)iterator.next();
            name = ((Number)feature.getAttribute(idField)).intValue();
            coordinate = ((Geometry)feature.getDefaultGeometry()).getCentroid().getCoordinate();
            double z = 0.0D;
            coordinate.z = z;
        }

        break MISSING_BLOCK_LABEL_117;
        Exception exception;
        exception;
        iterator.close();
        throw exception;
        iterator.close();
        return id2CoordinatesMap;
    }

    private void storeResult(double resultSwe[], double resultMel[], HashMap interpolatedCoordinatesMap)
        throws MismatchedDimensionException, Exception {
        WritableRandomIter outIterSWE = RandomIterFactory.createWritable(outWRSWE, null);
        WritableRandomIter outIterMEL = RandomIterFactory.createWritable(outWRMEL, null);
        Set pointsToInterpolateIdSett = interpolatedCoordinatesMap.keySet();
        Iterator idIterator = pointsToInterpolateIdSett.iterator();
        int c = 0;
        MathTransform transf = inRainGrid.getGridGeometry().getCRSToGrid2D();
        DirectPosition gridPoint = new DirectPosition2D();
        while(idIterator.hasNext())  {
            int id = ((Integer)idIterator.next()).intValue();
            Coordinate coordinate = (Coordinate)interpolatedCoordinatesMap.get(Integer.valueOf(id));
            DirectPosition point = new DirectPosition2D(inRainGrid.getCoordinateReferenceSystem(), coordinate.x, coordinate.y);
            transf.transform(point, gridPoint);
            double gridCoord[] = gridPoint.getCoordinate();
            int x = (int)gridCoord[0];
            int y = (int)gridCoord[1];
            outIterSWE.setSample(x, y, 0, resultSwe[c]);
            outIterMEL.setSample(x, y, 0, resultMel[c]);
            c++;
        }
        RegionMap regionMap = CoverageUtilities.gridGeometry2RegionParamsMap(inRainGrid.getGridGeometry());
        outMeltingDataGrid = CoverageUtilities.buildCoverage("gridded", outWRMEL, regionMap, inRainGrid.getGridGeometry().getCoordinateReferenceSystem());
        outSweDataGrid = CoverageUtilities.buildCoverage("gridded", outWRSWE, regionMap, inRainGrid.getGridGeometry().getCoordinateReferenceSystem());
    }

    private double getDeclination(double dayangb) {
        double delta = ((((0.37230000000000002D + 23.256699999999999D * Math.sin(dayangb)) - 0.75800000000000001D * Math.cos(dayangb)) + 0.1149D * Math.sin(2D * dayangb) + 0.36559999999999998D * Math.cos(2D * dayangb)) - 0.17119999999999999D * Math.sin(3D * dayangb)) + 0.0201D * Math.cos(3D * dayangb);
        return Math.toRadians(delta);
    }

    public GridCoverage2D inSkyview;
    public double pm1;
    public int pDimMeas;
    public String tStartDate;
    public int inTimestep;
    public String tEndDate;
    public GridCoverage2D inDem;
    public String time;
    public double pTmelt;
    public double pR;
    public GridCoverage2D inTempGrid;
    public GridCoverage2D incondIniIGrid;
    public GridCoverage2D incondIniLGrid;
    public GridCoverage2D inRainGrid;
    public GridCoverage2D inInsFeb;
    public GridCoverage2D inInsJan;
    public GridCoverage2D inInsMar;
    public GridCoverage2D inInsApr;
    public GridCoverage2D inInsMay;
    public GridCoverage2D inInsJun;
    public SimpleFeatureCollection inStations;
    public String fStationsid;
    public boolean doRaster;
    public HashMap inTemp;
    public HashMap inRainfall;
    public double pCmf;
    public double pCrf;
    public double pCff;
    public boolean doDaily;
    public double pCr;
    public double pCs;
    public String pathTemp;
    public int numSitesWhereCalibrate;
    public String pathToSolarRad;
    public String pathRainf;
    public int pMode;
    public boolean pDoReadMeas;
    public String pathRainfMaps;
    public String pathTempMaps;
    public IJGTProgressMonitor pm;
    public GridCoverage2D outGrid;
    public String pathToSwe;
    public String pathToPrainOutput;
    public String pathToSolidWater;
    public String pathToLiquidWater;
    public String pPathtoMeas;
    public double outMeasured[];
    public String pathToMelting;
    public GridCoverage2D outMeltingDataGrid;
    public GridCoverage2D outSweDataGrid;
    public double outSwevector[];
    private static final double pLapse = -0.0064999999999999997D;
    private HashMap attribute;
    private int height;
    private int width;
    public HashMap outMeltingData;
    public HashMap outSWEData;
    public HashMap outPrainData;
    private WritableRaster insFebWR;
    private WritableRaster insJanWR;
    private WritableRaster intempWR;
    private WritableRaster inrainWR;
    private WritableRaster demwr;
    private WritableRaster insMarWR;
    private WritableRaster insAprWR;
    private WritableRaster insJunWR;
    private WritableRaster insMayWR;
    public double lambda;
    private double lambdaVect[];
    private double xStation[];
    private double yStation[];
    private int idStation[];
    private int colnumvetVect[];
    private int rownumvetVect[];
    private boolean flag;
    private WritableRaster skyviewfactorWR;
    private HashMap rain_values;
    private HashMap solar_values;
    private HashMap temp_values;
    private HashMap condiniI_values;
    private HashMap condiniL_values;
    boolean init2;
    private double hoursToJan;
    private double hoursToFeb;
    private double hoursToMar;
    private double hoursToApr;
    private double hoursToMay;
    private double hoursToJune;
    public HashMap outMelting;
    public HashMap outSWE;
    private WritableRaster outWRSWE;
    private WritableRaster outWRMEL;
    private double outMELTINGVECTOR[];
    private double outSWEVECTOR[];
    private double outMELTINGVECTOR_forMaps[];
    private double outSWEVECTOR_forMaps[];
    public HashMap inCondIniI;
    public HashMap inCondIniL;
    public File folder;
    public File P_listOfFiles[];
    public File T_listOfFiles[];
    public int conta;
    private LinkedHashMap pointsToComputeId2Coordinates;
    private int cols;
    private int rows;
    private double south;
    private double west;
    private double xres;
    private double yres;
    private boolean doOneTime;
}
*/