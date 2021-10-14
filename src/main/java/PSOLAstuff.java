import com.github.psambit9791.jdsp.filter.Butterworth;
import com.github.psambit9791.jdsp.signal.peaks.FindPeak;
import com.github.psambit9791.jdsp.signal.peaks.Peak;
import com.github.psambit9791.jdsp.signal.peaks.Spike;
import com.github.psambit9791.jdsp.windows.Boxcar;
import com.github.psambit9791.jdsp.windows.Hanning;

import java.util.ArrayList;
import java.util.Arrays;

public class PSOLAstuff {
    public static double[] PSOLA(double[] signal, long sampleRate){
        FindPeak fp=new FindPeak(signal);
        Spike out=fp.getSpikes();
        Peak outpeak=fp.detectPeaks();
        int[] outMeanFilter=   outpeak.filterByProminence(20.0, 10000.0);
        //   int[] outMeanFilter = out.filterByProperty(1.0, 5000.0, "mean");
        //  int[] outMeanFilter = fp.detectRelativeMaxima();
        double[] spikeHeights=new double[outMeanFilter.length];
        double distsum=0;
        int co=0;
        double[] outMeanPoints=new double[outMeanFilter.length];
        double mindist=1000;
        for (int i =0; i<outMeanFilter.length; i++){

            if(signal[outMeanFilter[i]]>0){
                if(i!=0 && Math.abs(outMeanFilter[i]-outMeanFilter[i-1])<=mindist){

                    if(outMeanFilter[i-1]>=outMeanFilter[i]){continue;}
                    else{ outMeanFilter[i-1]=0;}
               //     System.out.println(outMeanFilter[i-1]);

                }
                outMeanPoints[i]=(double) outMeanFilter[i];
                spikeHeights[i]= signal[outMeanFilter[i]];

                if(i!=0&&outMeanPoints[i-1]!=0&&outMeanPoints[i]!=0){
               //     System.out.println("Hi");
                    for(int j=(i-1); j>=0; j--){
                        if(outMeanPoints[j]!=0){
                            distsum+=outMeanPoints[i]-outMeanPoints[j];
                            co++;
                            break;
                        }
                    }


                }}

        }
        double period=distsum/co;
        Boxcar hanning=new Boxcar((int) period, false);
       double[] window= hanning.getWindow();
        //double[] singlePeriod=getUnderAndModify(window, signal, 2000); -- 0 -- 0 -- 0 -- 0 -- 0 -- 0 -- 0 -- 0 -- 0 -- 0 -- 0 -- 0 -- 0 -- 0 -- 0 -- 0 -- 0 -- 0 -
        ArrayList<double[]> periods=getAllUnder(window, signal);
        ArrayList<double[]> shiftedPeriods=new ArrayList<>();
        for(int i=0; i<periods.size(); i++){
            shiftedPeriods.add(TribuoTest.shiftArray(periods.get(i), (int)(periods.get(i).length+(((period)/50)*(i)))));
            System.out.println((int) ((period/3)*(i+1)));
        }
        double[] sum=new double[signal.length];
        for (int i=0; i<shiftedPeriods.size(); i++){
            sum=addSignals(shiftedPeriods.get(i), sum);
        }
        System.out.println(period+"aa"+ window.length);
       for (int i=0; i< window.length; i++){
          window[i]=window[i]*2000;
       }
       double[] repeatedWindow=new double[signal.length];
       for (int i=0; i<signal.length; i++){
           repeatedWindow[i]=window[i% window.length];
       }

        //Butterworth bw=new Butterworth(sum, sampleRate);
       //bw.lowPassFilter(5, 20);
      //System.out.println(Arrays.toString(window));

       return sum;


    }
    public static double[] findSpikes(double[] signal, double distance){
    return null;
    }
    public static double[] getUnderAndModify(double[] window, double[] signal, int multiplier, int shift){
        double[] under=new double[signal.length];
        for (int i=shift; i< window.length+shift; i++) {
            under[i]=Math.min(signal[i], window[i-shift]*multiplier);

        }
        return under;
    }
    public static ArrayList<double[]> getAllUnder(double[] window, double[] signal){
        ArrayList<double[]> overlappedValues=new ArrayList<>();
        for (int i=0; i< Math.floor((signal.length)/(window.length)); i++){
            overlappedValues.add(getUnderAndModify(window, signal, 5000, i*window.length));
        }

        return overlappedValues;
    }

    public static double[] addSignals(double[] signal1, double[] signal2){
        double[] signalSum=new double[Math.min(signal1.length, signal2.length)];
    for (int i=0; i< Math.min(signal1.length, signal2.length); i++){
        signalSum[i]=signal1[i]+signal2[i];
    }
    return signalSum;
    }
}
