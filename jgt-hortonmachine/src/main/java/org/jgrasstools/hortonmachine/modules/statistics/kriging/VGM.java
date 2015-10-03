package org.jgrasstools.hortonmachine.modules.statistics.kriging;

import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;

public class VGM extends JGTModel {

    public VGM() {
        distances = null;
        Np = null;
        pm = new LogProgressMonitor();
        result = null;
        obs = null;
        inp = null;
    }

    public void process()
        throws Exception {
        result = calculate(distances, modelname, sill, range, nugget);
        obs = inp;
        Npairs = Np;
    }

    public static double[] calculate(double distance[], String model, double sill, double range, double nug) {
        double result[] = null;
        if("exponential".equals(model))
            result = fn_exponential(distance, sill, range, nug);
        if("gaussian".equals(model))
            result = fn_gaussian(distance, sill, range, nug);
        if("spherical".equals(model))
            result = fn_spherical(distance, sill, range, nug);
        if("pentaspherical".equals(model))
            result = fn_pentaspherical(distance, sill, range, nug);
        if("linear".equals(model))
            result = fn_linear(distance, sill, range, nug);
        if("circular".equals(model))
            result = fn_circolar(distance, sill, range, nug);
        if("bessel".equals(model))
            result = fn_bessel(distance, sill, range, nug);
        if("periodic".equals(model))
            result = fn_periodic(distance, sill, range, nug);
        if("hole".equals(model))
            result = fn_hole(distance, sill, range, nug);
        if("logaritmic".equals(model))
            result = fn_logaritmic(distance, sill, range, nug);
        if("power".equals(model))
            result = fn_power(distance, sill, range, nug);
        if("spline".equals(model))
            result = fn_spline(distance, sill, range, nug);
        return result;
    }

    double[] fn_nugget(double dist[], double range) {
        double nugget[] = new double[dist.length];
        for(int i = 0; i < dist.length; i++)
            nugget[i] = dist[i] != 0.0D ? 1.0D : 0.0D;

        return nugget;
    }

    public static double[] fn_exponential(double dist[], double sill, double range, double nug) {
        int length = dist.length;
        double func[] = new double[length];
        for(int i = 0; i < length; i++)
            if(dist[i] != 0.0D)
                func[i] = nug + sill * (1.0D - Math.exp(-dist[i] / range));

        return func;
    }

    public static double[] fn_gaussian(double dist[], double sill, double range, double nug) {
        int length = dist.length;
        double func[] = new double[length];
        for(int i = 0; i < length; i++) {
            double hr = dist[i] / range;
            if(dist[i] != 0.0D)
                func[i] = nug + sill * (1.0D - Math.exp(-(hr * hr)));
        }

        return func;
    }

    public static double[] fn_spherical(double dist[], double sill, double range, double nug) {
        int length = dist.length;
        double func[] = new double[length];
        for(int i = 0; i < length; i++) {
            double hr = dist[i] / range;
            if(dist[i] != 0.0D)
                func[i] = nug + sill * hr * (1.5D - 0.5D * hr * hr);
            if(dist[i] >= range)
                func[i] = sill;
        }

        return func;
    }

    public static double[] fn_pentaspherical(double dist[], double sill, double range, double nug) {
        int length = dist.length;
        double func[] = new double[length];
        double hr = 0.0D;
        for(int i = 0; i < length; i++) {
            hr = dist[i] / range;
            double h2r2 = hr * hr;
            if(dist[i] != 0.0D)
                func[i] = nug + sill * (hr * (1.875D + h2r2 * (-1.25D + h2r2 * 0.375D)));
            if(dist[i] >= range)
                func[i] = sill;
        }

        return func;
    }

    public static double[] fn_linear(double dist[], double sill, double range, double nug) {
        int length = dist.length;
        double func[] = new double[length];
        for(int i = 0; i < length; i++) {
            if(dist[i] != 0.0D)
                func[i] = nug + sill * (dist[i] / range);
            if(dist[i] >= range)
                func[i] = sill;
        }

        return func;
    }

    public static double[] fn_circolar(double dist[], double sill, double range, double nug) {
        int length = dist.length;
        double func[] = new double[length];
        for(int i = 0; i < length; i++) {
            double hr = dist[i] / range;
            if(dist[i] != 0.0D)
                func[i] = nug + sill * (0.63661977236758138D * (hr * Math.sqrt(1.0D - hr * hr) + Math.asin(hr)));
            if(dist[i] >= range)
                func[i] = sill;
        }

        return func;
    }

    public static double[] fn_bessel(double dist[], double sill, double range, double nug) {
        double MIN_BESS = 0.001D;
        double func[] = new double[dist.length];
        for(int i = 0; i < dist.length; i++) {
            double hr = dist[i] / range;
            if(hr > MIN_BESS)
                func[i] = nug + sill * (1.0D - hr * bessk1(hr));
        }

        return func;
    }

    static double bessk1(double x) {
        double ans;
        if(x <= 2D) {
            double y = (x * x) / 4D;
            ans = Math.log(x / 2D) * bessi1(x) + (1.0D / x) * (1.0D + y * (0.15443144D + y * (-0.67278579000000005D + y * (-0.18156897D + y * (-0.019194019999999999D + y * (-0.0011040399999999999D + y * -4.6860000000000002E-05D))))));
        } else {
            double y = 2D / x;
            ans = (Math.exp(-x) / Math.sqrt(x)) * (1.2533141400000001D + y * (0.23498619000000001D + y * (-0.036556199999999997D + y * (0.015042679999999999D + y * (-0.0078035300000000004D + y * (0.0032561399999999998D + y * -0.00068245000000000003D))))));
        }
        return ans;
    }

    static double bessi1(double x) {
        double ax;
        double ans;
        if((ax = Math.abs(x)) < 3.75D) {
            double y = x / 3.75D;
            y *= y;
            ans = ax * (0.5D + y * (0.87890594D + y * (0.51498869000000003D + y * (0.15084934D + y * (0.026587329999999999D + y * (0.0030153200000000002D + y * 0.00032411000000000001D))))));
        } else {
            double y = 3.75D / ax;
            ans = 0.02282967D + y * (-0.028953119999999999D + y * (0.01787654D - y * 0.0042005899999999997D));
            ans = 0.39894227999999998D + y * (-0.039880239999999997D + y * (-0.0036201800000000002D + y * (0.0016380100000000001D + y * (-0.01031555D + y * ans))));
            ans *= Math.exp(ax) / Math.sqrt(ax);
        }
        return x >= 0.0D ? ans : -ans;
    }

    public static double[] fn_periodic(double dist[], double sill, double range, double nug) {
        double func[] = new double[dist.length];
        for(int i = 0; i < dist.length; i++)
            if(dist[i] != 0.0D)
                func[i] = nug + sill * (1.0D - Math.cos((6.2831853071795862D * dist[i]) / range));

        return func;
    }

    public static double[] fn_hole(double dist[], double sill, double range, double nug) {
        double func[] = new double[dist.length];
        for(int i = 0; i < dist.length; i++)
            if(dist[i] != 0.0D)
                func[i] = nug + sill * (1.0D - Math.sin(dist[i] / range) / (dist[i] / range));

        return func;
    }

    public static double[] fn_logaritmic(double dist[], double sill, double range, double nug) {
        double func[] = new double[dist.length];
        for(int i = 0; i < dist.length; i++)
            if(dist[i] != 0.0D)
                func[i] = nug + sill * Math.log(dist[i] / range);

        return func;
    }

    public static double[] fn_power(double dist[], double sill, double range, double nug) {
        double func[] = new double[dist.length];
        for(int i = 0; i < dist.length; i++)
            if(dist[i] != 0.0D)
                func[i] = nug + sill * Math.pow(dist[i], range);

        return func;
    }

    public static double[] fn_spline(double dist[], double sill, double range, double nug) {
        double func[] = new double[dist.length];
        for(int i = 0; i < dist.length; i++) {
            if(dist[i] != 0.0D)
                func[i] = nug + sill * (dist[i] * dist[i] * Math.log(dist[i]));
            if(dist[i] >= range)
                func[i] = sill;
        }

        return func;
    }

    public double distances[];
    public double sill;
    public double Np[];
    public double Npairs[];
    public double range;
    public double nugget;
    public String modelname;
    public IJGTProgressMonitor pm;
    public double result[];
    public double obs[];
    public double inp[];
}
