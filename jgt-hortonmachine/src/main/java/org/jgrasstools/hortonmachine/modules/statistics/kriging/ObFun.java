package org.jgrasstools.hortonmachine.modules.statistics.kriging;


public class ObFun {

    public ObFun() {
        outObFunNash = 0.0D;
        outObFunRMSE = 0.0D;
        outObFunRMSEW = 0.0D;
        outObFunPBIAS = 0.0D;
        outR = 0.0D;
        outsdS_sdO = 0.0D;
        w = null;
        outmS_mO = 0.0D;
        outFHF = 0.0D;
        outFLF = 0.0D;
        outObFunFsls = 0.0D;
        outKGE = 0.0D;
    }

    public void process()
        throws Exception {
        double coef1_numO_S = 0.0D;
        double coef1_numO_S_W = 0.0D;
        double coef1_numS_O = 0.0D;
        double coef1_den = 0.0D;
        int contamedia = 0;
        double sommamediaoss = 0.0D;
        double sommamediasim = 0.0D;
        for(int i = 0; i < observed.length; i++)
            if(observed[i] != NaN && simulated[i] != NaN) {
                contamedia++;
                sommamediaoss += observed[i];
                sommamediasim += simulated[i];
            }

        double slsS_O = 0.0D;
        double slsO_S = 0.0D;
        double mediaoss = sommamediaoss / (double)contamedia;
        double mediasim = sommamediasim / (double)contamedia;
        int conta = 0;
        double numvaprev = 0.0D;
        double sommaosservati2 = 0.0D;
        double sumproduct = 0.0D;
        double sommasimulati2 = 0.0D;
        double sommaosservatimenoprevisti = 0.0D;
        double sommaprevistimenoosservati = 0.0D;
        double numR = 0.0D;
        double den1R = 0.0D;
        double den2R = 0.0D;
        double FLF = 0.0D;
        double FHF = 0.0D;
        double sommapesi = 0.0D;
        for(int i = 0; i < observed.length; i++)
            if(observed[i] != NaN && simulated[i] != NaN) {
                conta++;
                sumproduct += simulated[i] * observed[i];
                sommaosservati2 += observed[i] * observed[i];
                sommasimulati2 += simulated[i] * simulated[i];
                slsS_O += (simulated[i] - observed[i]) * (simulated[i] - observed[i]);
                slsO_S += (observed[i] - simulated[i]) * (observed[i] - simulated[i]);
                coef1_numO_S += (observed[i] - simulated[i]) * (observed[i] - simulated[i]);
                if(w != null) {
                    coef1_numO_S_W += (observed[i] - simulated[i]) * (observed[i] - simulated[i]) * w[i];
                    sommapesi += w[i];
                } else {
                    coef1_numO_S_W = coef1_numO_S + (observed[i] - simulated[i]) * (observed[i] - simulated[i]);
                    sommapesi++;
                }
                coef1_numS_O += (simulated[i] - observed[i]) * (simulated[i] - observed[i]);
                coef1_den += (observed[i] - mediaoss) * (observed[i] - mediaoss);
                numvaprev += (simulated[i] - mediasim) * (simulated[i] - mediasim);
                sommaosservatimenoprevisti += observed[i] - simulated[i];
                sommaprevistimenoosservati += simulated[i] - observed[i];
                numR += (observed[i] - mediaoss) * (simulated[i] - mediasim);
                den1R += (observed[i] - mediaoss) * (observed[i] - mediaoss);
                den2R += (simulated[i] - mediasim) * (simulated[i] - mediasim);
                if(simulated[i] != 0.0D && observed[i] != 0.0D)
                    FLF += (Math.log(simulated[i]) - Math.log(observed[i])) * (Math.log(simulated[i]) - Math.log(observed[i]));
                FHF += (simulated[i] - observed[i]) * (simulated[i] - observed[i]);
            }

        outFHF = FHF / (double)conta;
        outFLF = FLF / (double)conta;
        double varianzaosservati = coef1_den / (double)(conta - 1);
        double varianzasimulati = numvaprev / (double)(conta - 1);
        double sdsimulati = Math.sqrt(varianzasimulati);
        double sdosservati = Math.sqrt(varianzaosservati);
        double den11R = Math.sqrt(den1R);
        double den22R = Math.sqrt(den2R);
        double R = numR / (den11R * den22R);
        double alpha = sdsimulati / sdosservati;
        double beta = mediasim / mediaoss;
        double EDradicando = (R - 1.0D) * (R - 1.0D) + (alpha - 1.0D) * (alpha - 1.0D) + (beta - 1.0D) * (beta - 1.0D);
        outKGE = 1.0D - Math.sqrt(EDradicando);
        double Nash_SutcliffeS_O = 1.0D - coef1_numS_O / coef1_den;
        outObFunNash = Nash_SutcliffeS_O;
        double RMSEO_S = coef1_numO_S / (double)conta;
        double RMSEO_S_W = coef1_numO_S_W;
        outObFunRMSE = Math.sqrt(RMSEO_S);
        outObFunRMSEW = RMSEO_S_W;
        outObFunPBIAS = (sommaosservatimenoprevisti * 100D) / sommamediaoss;
        outObFunFsls = slsS_O / (double)conta;
        outR = R;
        outsdS_sdO = sdsimulati - sdosservati;
        outmS_mO = mediasim - mediaoss;
    }

    public double observed[];
    public double simulated[];
    public double NaN;
    public double outObFunNash;
    public double outObFunRMSE;
    public double outObFunRMSEW;
    public double outObFunPBIAS;
    public double outR;
    public double outsdS_sdO;
    public double w[];
    public double outmS_mO;
    public double outFHF;
    public double outFLF;
    public double outObFunFsls;
    public double outKGE;
}
