package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.hymod;

import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.IAdigeEngine;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.core.*;
import org.jgrasstools.hortonmachine.modules.network.PfafstetterNumber;
import org.joda.time.DateTime;
import org.opengis.feature.simple.SimpleFeature;

public class HymodAdigeEngineModified
    implements IAdigeEngine {

    public HymodAdigeEngineModified(Double pAlpha, Double pB, Double pCmax, boolean pDoGeom, boolean pDorouting, Double LAMBDA1, Double LAMBDA2, 
            Double QR, String fIdDistanze, Double vR, Double pRq, Double pRs, Double pMrain, Double pMetp, 
            String fIdnetnum, SimpleFeatureCollection inSubbasinDist, Double pQ0, List orderedHillslopes, HashMap index2Basinid, HashMap outDischarge, HashMap outSubDischarge, 
            HashMap outVelocity, List pfaffsList, boolean doLog, boolean doPrint, IJGTProgressMonitor pm, double initalconditios[]) {
        xSlow = null;
        xLoss = null;
        xQuick = null;
        dischargeContributorList = new ArrayList();
        conta = 0;
        ttTimestep = 0;
        this.pAlpha = pAlpha;
        this.pB = pB;
        this.pCmax = pCmax;
        this.pDoGeom = pDoGeom;
        this.pDorouting = pDorouting;
        this.LAMBDA1 = LAMBDA1;
        this.LAMBDA2 = LAMBDA2;
        this.QR = QR;
        this.fIdDistanze = fIdDistanze;
        this.vR = vR;
        this.pRq = pRq;
        this.pRs = pRs;
        this.pMrain = pMrain;
        this.pMetp = pMetp;
        this.fIdnetnum = fIdnetnum;
        this.inSubbasinDist = inSubbasinDist;
        this.pQ0 = pQ0;
        outActualEtp = outVelocity;
        this.orderedHillslopes = orderedHillslopes;
        this.index2Basinid = index2Basinid;
        this.outDischarge = outDischarge;
        this.outSubDischarge = outSubDischarge;
        this.pfaffsList = pfaffsList;
        this.doLog = doLog;
        this.doPrint = doPrint;
        this.pm = pm;
    }

    public void addDischargeContributor(IDischargeContributor dischargeContributor) {
        dischargeContributorList.add(dischargeContributor);
    }

    public HashMap getDischarge() {
        return outDischarge;
    }

    public HashMap getSubDischarge() {
        return outSubDischarge;
    }

    public HashMap getVelocity() {
        return outActualEtp;
    }

    public double[] solve(DateTime currentTimstamp, int tTimestep, double internalTimestepInMinutes, double initialConditions[], double rainArray[], double etpArray[])
        throws IOException {
        FeatureIterator iterator;
        ttTimestep = tTimestep;
        netnumList = new ArrayList();
        distanceList = new ArrayList();
       // if(!pDoGeom)
         //   break MISSING_BLOCK_LABEL_225;
        maxdist = 0.0D;
        mindist = 10000D;
      //  if(inSubbasinDist == null)
        //    break MISSING_BLOCK_LABEL_225;
        netnum2Dist = new HashMap();
        iterator = inSubbasinDist.features();
        while(iterator.hasNext())  {
            SimpleFeature feature = (SimpleFeature)iterator.next();
            int name = ((Number)feature.getAttribute(fIdnetnum)).intValue();
            double z = 0.0D;
            if(fIdDistanze != null)
                try {
                    z = ((Number)feature.getAttribute(fIdDistanze)).doubleValue();
                }
                catch(NullPointerException nullpointerexception) { }
            if(z > maxdist)
                maxdist = z;
            if(z < mindist)
                mindist = z;
            netnum2Dist.put(Integer.valueOf(name), Double.valueOf(z));
        }
       // break MISSING_BLOCK_LABEL_218;
     //   Exception exception;
    //    exception;
        iterator.close();
       //throw exception;
        iterator.close();
        if(initialConditions != null) {
            for(int i = orderedHillslopes.size() - 1; i >= 0; i--) {
                xLoss[i] = initialConditions[i];
                xSlow[i] = initialConditions[i + orderedHillslopes.size()];
                xQuick[0][i] = initialConditions[i + 2 * orderedHillslopes.size()];
                xQuick[1][i] = initialConditions[i + 3 * orderedHillslopes.size()];
                xQuick[2][i] = initialConditions[i + 4 * orderedHillslopes.size()];
            }

        }
        if(initialConditions == null) {
            xSlow = new double[orderedHillslopes.size()];
            coeffs = new double[orderedHillslopes.size()];
            xQuick = new double[3][orderedHillslopes.size()];
            xLoss = new double[orderedHillslopes.size()];
            initialConditions = new double[orderedHillslopes.size() * 6];
            for(int i = orderedHillslopes.size() - 1; i >= 0; i--) {
                IHillSlope hillSlope = (IHillSlope)orderedHillslopes.get(i);
                double areaKm2 = hillSlope.getHillslopeArea() / 1000000D;
                initialConditions[i + 5 * orderedHillslopes.size()] = pQ0.doubleValue() / (double)(i + 1);
                coeffs[i] = (Math.pow(10D, 9D) * (double)tTimestep * 60D) / (areaKm2 * Math.pow(10D, 12D));
                int idd = ((Integer)index2Basinid.get(Integer.valueOf(i))).intValue();
                double RRS = 0.0D;
                if(pDoGeom)
                    RRS = pRs.doubleValue() / (((Double)netnum2Dist.get(Integer.valueOf(idd))).doubleValue() / maxdist);
                else
                    RRS = pRs.doubleValue();
                xSlow[i] = (pQ0.doubleValue() * coeffs[i]) / RRS;
            }

        }
        for(int i = orderedHillslopes.size() - 1; i >= 0; i--) {
            IHillSlope hillSlope = (IHillSlope)orderedHillslopes.get(i);
            PfafstetterNumber pfaf = hillSlope.getPfafstetterNumber();
           // PfafstetterNumber pfaf = hillSlope.getPfafstetterNumber();
            double rain = rainArray[i];
            double etp = etpArray[i];
            if(rain < 0.0D)
                rain = 0.0D;
            if(etp < 0.0D)
                etp = 0.0D;
            rain *= pMrain.doubleValue();
            double out_excess[] = excess(xLoss[i], rain, etp);
            double UT1 = out_excess[0];
            double UT2 = out_excess[1];
            xLoss[i] = out_excess[2];
            Integer basinId = (Integer)index2Basinid.get(Integer.valueOf(i));
            double UQ = pAlpha.doubleValue() * UT2 + UT1;
            double US = (1.0D - pAlpha.doubleValue()) * UT2;
            int iddd = ((Integer)index2Basinid.get(Integer.valueOf(i))).intValue();
            double inflow = US;
            double RRRS = 0.0D;
            if(pDoGeom) {
                RRRS = pRs.doubleValue() / (((Double)netnum2Dist.get(Integer.valueOf(iddd))).doubleValue() / maxdist);
                if(RRRS > 0.10000000000000001D) {
                    RRRS = 0.10000000000000001D;
                    pRs = Double.valueOf(0.10000000000000001D * (((Double)netnum2Dist.get(Integer.valueOf(iddd))).doubleValue() / maxdist));
                }
                if(RRRS < 0.01D) {
                    RRRS = 0.01D;
                    pRs = Double.valueOf(0.01D * (((Double)netnum2Dist.get(Integer.valueOf(iddd))).doubleValue() / maxdist));
                }
            } else {
                RRRS = pRs.doubleValue();
            }
            double out_linres1[] = linres(xSlow[i], inflow, RRRS, 1.0D);
            xSlow[i] = out_linres1[0];
            double outflow1 = out_linres1[1];
            double QS = outflow1;
            inflow = UQ;
            double outflow2 = 0.0D;
            double RRRQ = 0.0D;
            if(pDoGeom) {
                RRRQ = pRq.doubleValue() / (((Double)netnum2Dist.get(Integer.valueOf(iddd))).doubleValue() / maxdist);
                if(RRRQ > 0.90000000000000002D) {
                    RRRQ = 0.90000000000000002D;
                    pRq = Double.valueOf(0.90000000000000002D * (((Double)netnum2Dist.get(Integer.valueOf(iddd))).doubleValue() / maxdist));
                }
                if(RRRQ < 0.10000000000000001D) {
                    RRRQ = 0.10000000000000001D;
                    pRq = Double.valueOf(0.10000000000000001D * (((Double)netnum2Dist.get(Integer.valueOf(iddd))).doubleValue() / maxdist));
                }
            } else {
                RRRQ = pRq.doubleValue();
            }
            for(int k = 0; k < 3; k++) {
                double out_linres2[] = linres(xQuick[k][i], inflow, RRRQ, 1.0D);
                xQuick[k][i] = out_linres2[0];
                outflow2 = out_linres2[1];
                inflow = outflow2;
            }

            double basinDischarge = (QS + outflow2) / coeffs[i];
            double basinSubDischarge = xLoss[i];
            double allContributionsDischarge = handleContributors(hillSlope, basinDischarge, i, initialConditions);
            basinDischarge = allContributionsDischarge;
            if(pfaffsList.contains(pfaf.toString())) {
                outDischarge.put(basinId, new double[] {
                    basinDischarge
                });
                outSubDischarge.put(basinId, new double[] {
                    basinSubDischarge
                });
                outActualEtp.put(basinId, new double[] {
                    etp
                });
            }
            initialConditions[i] = xLoss[i];
            initialConditions[i + orderedHillslopes.size()] = xSlow[i];
            initialConditions[i + 2 * orderedHillslopes.size()] = xQuick[0][i];
            initialConditions[i + 3 * orderedHillslopes.size()] = xQuick[1][i];
            initialConditions[i + 4 * orderedHillslopes.size()] = xQuick[2][i];
            initialConditions[i + 5 * orderedHillslopes.size()] = basinDischarge;
            outDischargeInternal.put(basinId, new double[] {
                basinDischarge
            });
        }

        return initialConditions;
    }

    private double handleContributors(IHillSlope hillSlope, double basinDischarge, int indice, double initialcond[]) {
        double summedContributions = 0.0D;
        List connectedUpstreamHillSlopes = hillSlope.getConnectedUpstreamElements();
        hillSlope.getConnectedDownstreamElement();
        double areamonte = ((IHillSlope)orderedHillslopes.get(indice)).getUpstreamArea(null);
        if(connectedUpstreamHillSlopes != null) {
            for(Iterator iterator = connectedUpstreamHillSlopes.iterator(); iterator.hasNext();) {
                IHillSlope tmpHillSlope = (IHillSlope)iterator.next();
                PfafstetterNumber pNum = tmpHillSlope.getPfafstetterNumber();
                int hillslopeId = tmpHillSlope.getHillslopeId();
                double upstreamDischarge = ((double[])outDischargeInternal.get(Integer.valueOf(hillslopeId)))[0];
                for(Iterator iterator1 = dischargeContributorList.iterator(); iterator1.hasNext();) {
                    IDischargeContributor dContributor = (IDischargeContributor)iterator1.next();
                    Double contributedDischarge = dContributor.getDischarge(pNum.toString());
                    if(!JGTConstants.isNovalue(contributedDischarge.doubleValue())) {
                        if(doLog && doPrint)
                            pm.message((new StringBuilder("----> For hillslope ")).append(hillSlope.getPfafstetterNumber()).append(" using hydrometer/dams data in pfafstetter: ").append(pNum.toString()).append("(meaning added ").append(contributedDischarge).append(" instead of ").append(upstreamDischarge).append(")").toString());
                        if(!JGTConstants.isNovalue(contributedDischarge.doubleValue()))
                            upstreamDischarge = dContributor.mergeWithDischarge(contributedDischarge.doubleValue(), upstreamDischarge);
                    }
                }

                summedContributions += upstreamDischarge;
            }

        }
        double output = 0.0D;
        boolean pdoRouting = pDorouting;
        if(pdoRouting)
            output = doRouting(summedContributions, basinDischarge, hillSlope, initialcond[indice + 5 * orderedHillslopes.size()], areamonte);
        else
            output = basinDischarge + summedContributions;
        return output;
    }

    private double doRouting(double UPdischarge, double HYMODbasinDischarge, IHillSlope hillslope, double condini, 
            double Areaamonte) {
        double y = 0.0D;
        int a = 0;
        int b = ttTimestep * 60;
        y = condini;
        double lambda1 = 0.0D;
        double lambda2 = 0.0D;
        double Qr = 0.0D;
        double vr = 0.0D;
        if(!JGTConstants.isNovalue(LAMBDA1.doubleValue()))
            lambda1 = LAMBDA1.doubleValue();
        else
            lambda1 = 0.20000000000000001D;
        if(!JGTConstants.isNovalue(LAMBDA2.doubleValue()))
            lambda2 = LAMBDA2.doubleValue();
        else
            lambda2 = -0.10000000000000001D;
        if(!JGTConstants.isNovalue(QR.doubleValue()))
            Qr = QR.doubleValue();
        else
            Qr = 100D;
        if(!JGTConstants.isNovalue(vR.doubleValue()))
            vr = vR.doubleValue();
        else
            vr = 1.5D;
        double Ar = 1.0D;
        double l = hillslope.getLinkLength();
        double Q_nMin1 = condini;
        int nmax = 1000;
        double t[] = new double[nmax];
        t[0] = 0.0D;
        double dt = 360D;
        int tfinal = ttTimestep * 60;
        for(int n = 0; n < nmax; n++) {
            if(dt + t[n] > (double)tfinal)
                dt = (double)tfinal - t[n];
            if(t[n] >= (double)tfinal)
                break;
            double result[] = newtonRouting(Q_nMin1, t[n], dt, UPdischarge, HYMODbasinDischarge, hillslope, lambda1, lambda2, vr, Qr, Ar, l);
            y = result[1];
            Q_nMin1 = y;
            t[n + 1] = t[n] + result[0];
        }

        vel = vr * Math.pow(y, lambda1) * Math.pow(hillslope.getHillslopeArea(), lambda2);
        return y;
    }

    private double computegRouting(double Q_n, double Q_nMin1, double t, double dt, double UPdischarge, double HYMODbasinDischarge, IHillSlope hillslope, double lambda1, double lambda2, double vr, double Qr, 
            double Ar, double l) {
        double yy = 0.0D;
        double func = EFFE(hillslope, lambda1, lambda2, vr, Qr, Ar, l, Q_n, HYMODbasinDischarge, UPdischarge);
        yy = (Q_nMin1 - Q_n) + dt * func;
        return yy;
    }

    private double dgRouting(double Q_n, double Q_nMin1, double t, double dt, double UPdischarge, double HYMODbasinDischarge, IHillSlope hillslope, double lambda1, double lambda2, double vr, double Qr, 
            double Ar, double l) {
        double der = 0.0D;
        double eps = 9.9999999999999995E-08D;
        der = (computegRouting(Q_n + eps, Q_nMin1, t, dt, UPdischarge, HYMODbasinDischarge, hillslope, lambda1, lambda2, vr, Qr, Ar, l) - computegRouting(Q_n - eps, Q_nMin1, t, dt, UPdischarge, HYMODbasinDischarge, hillslope, lambda1, lambda2, vr, Qr, Ar, l)) / (2D * eps);
        return der;
    }

    private double[] newtonRouting(double Q_nMin1, double t, double dt, double UPdischarge, double HYMODbasinDischarge, IHillSlope hillslope, double lambda1, double lambda2, double vr, double Qr, double Ar, 
            double l) {
        double r[] = new double[2];
        double Q_n = Q_nMin1;
        double tol = 9.9999999999999995E-08D;
        double nmax = 100D;
        double dtout = dt;
        double flag = 0.0D;
        int k = 0;
        for(k = 0; (double)k <= nmax; k++) {
            double ter1 = computegRouting(Q_n, Q_nMin1, t, dt, UPdischarge, HYMODbasinDischarge, hillslope, lambda1, lambda2, vr, Qr, Ar, l);
            if(Math.abs(ter1) < tol)
                break;
            double ter2 = dgRouting(Q_n, Q_nMin1, t, dt, UPdischarge, HYMODbasinDischarge, hillslope, lambda1, lambda2, vr, Qr, Ar, l);
            double dq = ter1 / ter2;
            double delta = 1.0D;
            int j = 0;
            for(j = 0; (double)j <= nmax; j++) {
                double ter3 = computegRouting(Q_n - dq * delta, Q_nMin1, t, dt, UPdischarge, HYMODbasinDischarge, hillslope, lambda1, lambda2, vr, Qr, Ar, l);
                double ter4 = computegRouting(Q_n, Q_nMin1, t, dt, UPdischarge, HYMODbasinDischarge, hillslope, lambda1, lambda2, vr, Qr, Ar, l);
                if(Math.abs(ter3) <= Math.abs(ter4))
                    break;
                delta /= 2D;
            }

            if((double)j >= nmax)
                System.out.println("Newton inner loop warning");
            Q_n -= delta * dq;
        }

        if(JGTConstants.isNovalue(Q_n)) {
            System.out.println("Newton produced NaN");
            flag = 1.0D;
        }
        if((double)k >= nmax) {
            System.out.println("Newton outer loop warning");
            flag = 1.0D;
        }
        if(flag == 1.0D) {
            System.out.println("Enter in bisection");
            double Qa = 0.0D;
            double Qb = 10D * Q_nMin1;
            for(int j = 0; (double)j <= nmax; j++) {
                double ter5 = computegRouting(Qa, Q_nMin1, t, dt, UPdischarge, HYMODbasinDischarge, hillslope, lambda1, lambda2, vr, Qr, Ar, l);
                double ter6 = computegRouting(Qb, Q_nMin1, t, dt, UPdischarge, HYMODbasinDischarge, hillslope, lambda1, lambda2, vr, Qr, Ar, l);
                if(ter5 * ter6 <= 0.0D)
                    break;
                System.out.println("Severe Bisection Error!!Initial Sign check falied");
                dt *= 0.10000000000000001D;
            }

            dtout = dt;
            for(int kk = 0; (double)kk < 10D * nmax; kk++) {
                double Qc = 0.5D * (Qa + Qb);
                double ter7 = computegRouting(Qa, Q_nMin1, t, dt, UPdischarge, HYMODbasinDischarge, hillslope, lambda1, lambda2, vr, Qr, Ar, l);
                double ter8 = computegRouting(Qc, Q_nMin1, t, dt, UPdischarge, HYMODbasinDischarge, hillslope, lambda1, lambda2, vr, Qr, Ar, l);
                if(ter7 * ter8 > 0.0D)
                    Qa = Qc;
                else
                    Qb = Qc;
                double res = Math.abs(computegRouting(Qc, Q_nMin1, ter8, dt, UPdischarge, HYMODbasinDischarge, hillslope, lambda1, lambda2, vr, Qr, Ar, l));
                if(res >= tol)
                    continue;
                System.out.println("Bisection OK!!!!!");
                Q_n = Qc;
                break;
            }

        }
        r[0] = dtout;
        r[1] = Q_n;
        return r;
    }

    private double T(IHillSlope hillslope, double areaamonte, double lambda1, double lambda2, 
            double vr, double Qr, double Ar, double l, double portata) {
        double re = 0.0D;
        double vvv = vr * Math.pow(portata, lambda1) * Math.pow(hillslope.getHillslopeArea(), lambda2);
        double kkk = vvv / (l * (1.0D - lambda1));
        if(kkk == 0.0D)
            re = 1.0000000000000001E-05D;
        else
            re = kkk;
        return re;
    }

    private double EFFE(IHillSlope hillslope, double lambda1, double lambda2, double vr, 
            double Qr, double Ar, double l, double portata, double HYMOD, double QTRIB) {
        double re = 0.0D;
        double uno = (QTRIB + HYMOD) - portata;
        double due = (vr * Math.pow(hillslope.getHillslopeArea(), lambda2) * Math.pow(portata, lambda1)) / (l * (1.0D - lambda1));
        if(due == 0.0D)
            due = 0.0001D;
        re = uno * due;
        return re;
    }

    private double[] excess(double x_losss, double Pval, double PETval) {
        double o_exces[] = new double[3];
        double xn_prev = x_losss;
        double coeff1 = 1.0D - ((pB.doubleValue() + 1.0D) * xn_prev) / pCmax.doubleValue();
        if(coeff1 < 0.0D) {
            System.out.println((new StringBuilder("coeff1= ")).append(coeff1).append("is assumed equal to 0").toString());
            coeff1 = 0.0D;
        }
        double exp = 1.0D / (pB.doubleValue() + 1.0D);
        double ct_prev = pCmax.doubleValue() * (1.0D - Math.pow(coeff1, exp));
        double UT1 = Math.max((Pval - pCmax.doubleValue()) + ct_prev, 0.0D);
        Pval -= UT1;
        if(Pval != Pval)
            System.out.println("ATTENTION: NaN for some system variable");
        double dummy = Math.min((ct_prev + Pval) / pCmax.doubleValue(), 1.0D);
        double coeff2 = 1.0D - dummy;
        double exp2 = pB.doubleValue() + 1.0D;
        double xn = (pCmax.doubleValue() / (pB.doubleValue() + 1.0D)) * (1.0D - Math.pow(coeff2, exp2));
        double UT2 = Math.max(Pval - (xn - xn_prev), 0.0D);
        double evap = Math.min(xn, PETval);
        double max = pCmax.doubleValue() / (pB.doubleValue() + 1.0D);
        evap = (1.0D - (max - xn) / max) * PETval;
        xn -= evap;
        PETval = evap;
        o_exces[0] = UT1;
        o_exces[1] = UT2;
        o_exces[2] = xn;
        if(xn != xn || UT1 != UT1 || UT2 != UT2)
            System.out.println("ATTENTION: NaN for some system variable");
        return o_exces;
    }

    private double[] linres(double x_sloww, double infloww, double RR, double dt) {
        double o_linres[] = new double[2];
        double x_sloww_prev = x_sloww;
        double x_slow_new = infloww * (1.0D - RR) + x_sloww_prev * (1.0D - RR);
        double outfloww = (x_slow_new * RR) / (1.0D - RR);
        o_linres[0] = x_slow_new;
        o_linres[1] = outfloww;
        return o_linres;
    }

    private final List orderedHillslopes;
    private double xSlow[];
    private double xLoss[];
    private double xQuick[][];
    private final HashMap outDischarge;
    private final HashMap outActualEtp;
    private final HashMap outSubDischarge;
    private final HashMap outDischargeInternal = new HashMap();
    private final HashMap index2Basinid;
    private List dischargeContributorList;
    private final boolean doPrint;
    private final boolean doLog;
    private final IJGTProgressMonitor pm;
    private final List pfaffsList;
    private double coeffs[];
    private List netnumList;
    private List distanceList;
    private HashMap netnum2Dist;
    private Double pAlpha;
    private Double pB;
    private Double pCmax;
    private boolean pDoGeom;
    private boolean pDorouting;
    private Double LAMBDA1;
    private Double LAMBDA2;
    private Double QR;
    private String fIdDistanze;
    private Double vR;
    private Double pRq;
    private Double pRs;
    private Double pMrain;
    private Double pMetp;
    private String fIdnetnum;
    private SimpleFeatureCollection inSubbasinDist;
    private Double pQ0;
    int conta;
    private double maxdist;
    private double mindist;
    private int ttTimestep;
    private double vel;
}
