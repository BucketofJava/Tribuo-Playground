import be.tarsos.dsp.AudioDispatcher;
import be.tarsos.dsp.AudioEvent;
import be.tarsos.dsp.PitchShifter;
import be.tarsos.dsp.WaveformSimilarityBasedOverlapAdd;
import be.tarsos.dsp.io.PipedAudioStream;
import be.tarsos.dsp.io.TarsosDSPAudioFormat;
import be.tarsos.dsp.io.TarsosDSPAudioInputStream;
import be.tarsos.dsp.io.jvm.AudioDispatcherFactory;
import be.tarsos.dsp.io.jvm.AudioPlayer;
import be.tarsos.dsp.io.jvm.JVMAudioInputStream;
import be.tarsos.dsp.writer.WriterProcessor;
import com.github.psambit9791.jdsp.filter.Butterworth;
import com.github.psambit9791.jdsp.io.Wav;
import com.github.psambit9791.jdsp.misc.Plotting;
import com.github.psambit9791.jdsp.signal.Detrend;
import com.github.psambit9791.jdsp.signal.Resample;
import com.github.psambit9791.jdsp.signal.peaks.FindPeak;
import com.github.psambit9791.jdsp.signal.peaks.Peak;
import com.github.psambit9791.jdsp.signal.peaks.Spike;
import com.github.psambit9791.jdsp.transform.DiscreteFourier;
import com.github.psambit9791.jdsp.transform.Hilbert;
import com.github.psambit9791.jdsp.transform.InverseDiscreteFourier;
import com.github.psambit9791.jdsp.windows.Hanning;
import com.github.psambit9791.wavfile.WavFileException;
import org.apache.commons.math3.complex.Complex;
import org.tribuo.*;
import org.tribuo.data.csv.CSVLoader;
import org.tribuo.datasource.ListDataSource;
import org.tribuo.evaluation.TrainTestSplitter;
import org.tribuo.math.optimisers.AdaGrad;
import org.tribuo.math.optimisers.SGD;
import org.tribuo.regression.*;
import org.tribuo.regression.evaluation.RegressionEvaluation;
import org.tribuo.regression.evaluation.RegressionEvaluator;
import org.tribuo.regression.rtree.CARTRegressionTrainer;
import org.tribuo.regression.sgd.linear.LinearSGDTrainer;
import org.tribuo.regression.sgd.objectives.SquaredLoss;
import org.tribuo.regression.xgboost.XGBoostRegressionTrainer;

import javax.sound.sampled.*;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;

public class TribuoTest {

    public static void main(String[] args) throws IOException, WavFileException, UnsupportedAudioFileException, LineUnavailableException {
        tarsosStuff("src/main/resources/Hello-SoundBible.com-218208532.wav");
        runWavStuff("src/main/resources/Hello-SoundBible.com-218208532.wav");
    }
    public static Model<Regressor> train(String name, Trainer<Regressor> trainer, Dataset<Regressor> trainData){
    Long startTime=System.currentTimeMillis();
    Model<Regressor> model=trainer.train(trainData);
    Long endTime=System.currentTimeMillis();
    System.out.println("Training took: "+(endTime-startTime)+"ms");
    RegressionEvaluator evaluator=new RegressionEvaluator();
    RegressionEvaluation evaluation=evaluator.evaluate(model, trainData);
    Regressor dimension=new Regressor("DIM-0", Double.NaN);
    System.out.printf("Evaluation (train):%n RMSE %f%n MAE %f%n R^2 %f%n", evaluation.rmse(dimension), evaluation.mae(dimension), evaluation.r2(dimension));
    return model;

    }
    public static void tarsosStuff(String wavname) throws IOException, WavFileException, UnsupportedAudioFileException, LineUnavailableException {

        Wav wav=new Wav();
        wav.readWav(wavname);
        Plotting fig=new Plotting();

        double[][] doubles=wav.getData("double");
        float[] buffer=new float[doubles.length];
        for (int i=0; i< doubles.length; i++){
            buffer[i]=(float) doubles[i][0];
        }
        final byte[] byteBuffer = new byte[buffer.length * 2];
        int bufferIndex = 0;
        for (int i = 0; i < byteBuffer.length; i++) {
            final int x = (int) (buffer[bufferIndex++] * 32767.0);
            byteBuffer[i] = (byte) x;
            i++;
            byteBuffer[i] = (byte) (x >>> 8);
        }
                AudioDispatcher adp =  AudioDispatcherFactory.fromFile(new File(wavname), byteBuffer.length, buffer.length);
        TarsosDSPAudioFormat formati=adp.getFormat();
        PitchShifter ps=new PitchShifter(1.1, formati.getSampleRate(), byteBuffer.length, buffer.length);
        adp.addAudioProcessor(ps);
         adp.addAudioProcessor(new AudioPlayer(JVMAudioInputStream.toAudioFormat(formati)));
          new Thread(adp).start();
//        Wav wav2=new Wav();
//        wav2.readWav("Hiyay.wav");
//        double[][] data2=wav2.getData("long");
//        double[] signal=new double[data2.length];
//        double[] indarr=new double[data2.length];
//        for(int i=0; i< data2.length; i++){
//            signal[i]=data2[i][0];
//            indarr[i]=i;
//        }
//        fig.addSignal("Cool", indarr, signal, false);
//        fig.plot();

    }
    public static void runWavStuff(String wavname) throws IOException, WavFileException, UnsupportedAudioFileException, LineUnavailableException {
        Wav wav=new Wav();
        wav.readWav(wavname);
        double[][] data=wav.getData("long");
        //System.out.println(Arrays.deepToString(data));
        int width = 600;
        int height = 500;
        String title = "Sample Figure";
        String x_axis = "Time";
        String y_axis = "Signal";
        Plotting fig = new Plotting(width, height, title, x_axis, y_axis);
        fig.initialisePlot();
        double[] xvals=new double[10000];
        double[] yvals=new double[10000];
        double[] carr=new double[10000];
        int c=0;
        double scalefactor=0.1;
        for(int i=4000; i<14000; i++){
            System.out.println(Arrays.toString(data[i]));
          //  if(Arrays.asList(xvals).contains(data[i][0])){ c++; continue;}

            carr[c]=c;
            xvals[c]=data[i][0]*scalefactor;
            yvals[c]=data[i][1]*scalefactor;

            c++;

        }
       // System.out.println(Arrays.toString(xvals));

        // For vertical lines
        double[][] verLines = {{0.0, 2.0, 3.4},
                {1.0, 4.0, 6.7},
                {2.0, 2.0, 2.2},
                {3.0, 3.0, 1.6},
                {4.0, 1.0, 3.6}};

      /*  for (int i=0; i<=verLines.length-1; i++) {
            fig.vline(verLines[i][0], verLines[i][1], verLines[i][2]);
        }
        */
// For horizontal lines
        fig.hline(0.0, 4.0, 4.0);
        fig.hline(0.0, 4.0, 6.7);

        Butterworth testButter=new Butterworth(yvals, (Long) wav.getProperties().get("SampleRate"));
        Butterworth testButter1=new Butterworth(xvals, (Long) wav.getProperties().get("SampleRate"));
        Detrend d1 = new Detrend(yvals, 3);
        double[] deterend = d1.detrendSignal();
        fig.addSignal("Signal", carr, yvals, false);
        FindPeak fp=new FindPeak(yvals);
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

            if(yvals[outMeanFilter[i]]>0){
                if(i!=0 && Math.abs(outMeanFilter[i]-outMeanFilter[i-1])<=mindist){

                    if(outMeanFilter[i-1]>=outMeanFilter[i]){continue;}
                    else{ outMeanFilter[i-1]=0;}
                    System.out.println(outMeanFilter[i-1]);

                }
            outMeanPoints[i]=(double) outMeanFilter[i];
            spikeHeights[i]= yvals[outMeanFilter[i]];

            if(i!=0&&outMeanPoints[i-1]!=0&&outMeanPoints[i]!=0){
                System.out.println("Hi");
                for(int j=(i-1); j>=0; j--){
                    if(outMeanPoints[j]!=0){
                        distsum+=outMeanPoints[i]-outMeanPoints[j];
                        co++;
                        break;
                    }
                }

            }}

        }
        System.out.println("Average is " + (distsum
                /co));

       double[] window= PSOLAstuff.PSOLA(yvals,(Long) wav.getProperties().get("SampleRate"));
        // fig.addPoints("Peaks", outMeanPoints, spikeHeights);
        Hanning hanning=new Hanning((int) carr.length, false);
        double[] window2= hanning.getWindow();
        String outputFileName = "signal.png";
       fig.addSignal("S", window, false);
// To plot on a window

// To save as an image
       // fig.saveAsPNG(outputFileName);
        Wav writeObj=new Wav();
        double[] lpx=testButter1.lowPassFilter(4, 1000);
        double[] lpy=testButter.lowPassFilter(4, 1000);
      //  System.out.println(yvals.length);
       // System.out.println(yvals.length*yvals.length);
        double[][] doubles=wav.getData("double");
        float[] buffer=new float[doubles.length];
        for (int i=0; i< doubles.length; i++){
            buffer[i]=(float) doubles[i][0];
        }
        final byte[] byteBuffer = new byte[buffer.length * 2];
        int bufferIndex = 0;
        for (int i = 0; i < byteBuffer.length; i++) {
            final int x = (int) (buffer[bufferIndex++] * 32767.0);
            byteBuffer[i] = (byte) x;
            i++;
            byteBuffer[i] = (byte) (x >>> 8);
        }
        File outFile = new File("out.wav");
        boolean bigEndian = false;
        boolean signed = true;
        int bits = 16;
        int channels = 1;
        //AudioFormat format;
        //format = new AudioFormat((float) (Long) wav.getProperties().get("SampleRate"), ((Long) wav.getProperties().get("ValidBits")).intValue(), ((Long) wav.getProperties().get("Channels")).intValue(), signed, bigEndian);
        // ByteArrayInputStream bais = new ByteArrayInputStream(byteBuffer);
       // AudioInputStream audioInputStream;
       // audioInputStream = new AudioInputStream(bais, format,buffer.length);
      //  AudioSystem.write(audioInputStream, AudioFileFormat.Type.WAVE, outFile);
      //  audioInputStream.close();
        //PipedAudioStream pas=new PipedAudioStream(wavname);
        //TarsosDSPAudioInputStream tais= pas.getMonoStream(((Long) wav.getProperties().get("SampleRate")).intValue(), 0);
        //byte[] bufferi=new byte[buffer.length*2];
        //tais.read(bufferi, 0, bufferi.length);
//        AudioDispatcher adp =  AudioDispatcherFactory.fromFile(new File(wavname), byteBuffer.length, buffer.length);
//        TarsosDSPAudioFormat formati=adp.getFormat();
//        PitchShifter ps=new PitchShifter(0.8, formati.getSampleRate(), byteBuffer.length, buffer.length);
//        adp.addAudioProcessor(ps);
        // adp.addAudioProcessor(new AudioPlayer(JVMAudioInputStream.toAudioFormat(formati)));
      //  new Thread(adp).start();

   /*  System.out.println("Starting fourier transform");
        DiscreteFourier df=new DiscreteFourier(yvals);
        System.out.println("Created Fourier object");
        df.dft(
        System.out.println("Ran dft");
        double[][] complexes=df.returnFull(false);
        System.out.println("Looping through values");
        double[] yyvals=new double[yvals.length];
        c=0;
        for(double[] complex:complexes){
            int cc=0;
            for(double d:complex){
            complex[cc]=d;
            cc++;}
            yyvals[c]=complex[1];
            c++;
        }

        System.out.println("Starting inverse Fourier");
       InverseDiscreteFourier idf=new InverseDiscreteFourier(shiftArray(complexes, -15000), false);
       InverseDiscreteFourier idf2=new InverseDiscreteFourier(shiftArray(complexes, 15000), false);
        idf.idft();
        idf2.idft();
   //   fig.addSignal("Signal 2", carr, yyvals, false);
       // fig.addSignal("Signal 3", carr, df.returnAbsolute(false), false);
        fig.addSignal("Signal 6", carr, idf.getRealSignal(), false);
     //   fig.addSignal("Signal 7", carr, idf2.getRealSignal(), false);
        System.out.println("Successfully added signal");
        DiscreteFourier df2=new DiscreteFourier(shiftArray(df.returnAbsolute(false), 5000));
        df2.dft();
//      fig.addSignal("Signal 5", carr, shiftArray(df2.returnAbsolute(false), 5000), false);
*/
      //  Hilbert h=new Hilbert(yvals);
        //h.hilbertTransform();
        //double[] frequency=Arrays.copyOf(h.getInstantaneousFrequency((Long) wav.getProperties().get("SampleRate")), carr.length);
       // System.out.println(carr.length);
        //System.out.println(frequency.length);
      //  fig.addSignal("Signal 4", carr, frequency, false);
        System.out.println(Arrays.deepToString(shiftArray(new double[][]{{1, 0},{2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}}, 3)));
      //  double[] fShiftedxvals=fourierShift(xvals, -100);
        //double[] fShiftedyvals=fourierShift(yvals, 0 );
        double[] fShiftedyvals=window;
        double[] fShiftedxvals=new double[fShiftedyvals.length];
       // fig.addSignal("Signal 8", carr, fShiftedxvals, false);
       //fig.addSignal("Signal 9", carr, fShiftedyvals, false);
     //   fig.addSignal("Signal 10", carr, lpy, false);

        double[][] newvals=new double[data.length][2];
        for (int i=0; i< xvals.length; i++){
            newvals[i]=new double[]{fShiftedxvals[i], fShiftedyvals[i]};
        }
        fig.plot();

        writeObj.putData(newvals, (Long) wav.getProperties().get("SampleRate"), "long", "Fun.wav");
        System.out.println("Data write successful");

    }

    public static double[] shiftArray(double[] doubles, int amount){
        double[] filler=new double[doubles.length];
        double[] returnval=new double[doubles.length];
        System.arraycopy(doubles, 0, filler, 0, filler.length);
        for(int i=0; i<doubles.length; i++){
            returnval[i]=filler[(i+amount)%filler.length];
        }
        return returnval;
    }
    public static double[][] shiftArray(double[][] doubles, int amount){
        double[][] filler=new double[doubles.length][doubles[0].length];
        double[][] returnval=new double[doubles.length][doubles[0].length];
        System.arraycopy(doubles, 0, filler, 0, filler.length);
        for(int i=0; i<doubles.length; i++){
            if(i+amount>=0){
            returnval[i]=filler[(i+amount)%filler.length];}
            else{
                returnval[i]=filler[((i+amount)%filler.length)+filler.length - 1];
            }
        }
        return returnval;
    }

    public static double[] fourierShift(double[] signal, int amount){
        System.out.println("Starting fourier transform");
        DiscreteFourier df=new DiscreteFourier(signal);
        System.out.println("Created Fourier object");
        df.dft();
        System.out.println("Ran dft");
        double[][] complexes=df.returnFull(false);
        System.out.println("Looping through values");
        double[] yyvals=new double[signal.length];
       int c=0;
        for(double[] complex:complexes){
            int cc=0;
            for(double d:complex){
                complex[cc]=d;
                cc++;}
            yyvals[c]=complex[1];
            c++;
        }

        System.out.println("Starting inverse Fourier");
       InverseDiscreteFourier idf=new InverseDiscreteFourier(shiftArray(complexes, amount), false);
       System.out.println("Inverse fourier initialized");
     //   InverseDiscreteFourier idf2=new InverseDiscreteFourier(shiftArray(complexes, 15000), false);
        idf.idft();
       // idf2.idft();
        //   fig.addSignal("Signal 2", carr, yyvals, false);
        // fig.addSignal("Signal 3", carr, df.returnAbsolute(false), false);

        //   fig.addSignal("Signal 7", carr, idf2.getRealSignal(), false);
        System.out.println("Successfully underwent inverse fourier");
       // DiscreteFourier df2=new DiscreteFourier(shiftArray(df.returnAbsolute(false), 5000));
       // df2.dft();
//      fig.addSignal("Signal 5", carr, shiftArray(df2.returnAbsolute(false), 5000), false);
        return idf.getRealSignal();
    }

    public static void evaluate(Model<Regressor> model, Dataset<Regressor> testData){
        RegressionEvaluator evaluator=new RegressionEvaluator();
        RegressionEvaluation evaluation=evaluator.evaluate(model, testData);
        Regressor dimension=new Regressor("DIM-0", Double.NaN);
        System.out.printf("Evaluation (test):%n RMSE %f%n MAE %f%n R^2 %f%n", evaluation.rmse(dimension), evaluation.mae(dimension), evaluation.r2(dimension));

    }
    public static void run() throws IOException {RegressionFactory regressionFactory=new RegressionFactory();
        CSVLoader<Regressor> loader=new CSVLoader<>(';', regressionFactory);
        ListDataSource<Regressor> source=loader.loadDataSource(Paths.get("src\\main\\resources\\winequality-red.csv"), "quality");
        TrainTestSplitter<Regressor> splitter=new TrainTestSplitter<>(source, 0.7f,0L);
        Dataset<Regressor> trainData= new MutableDataset<>(splitter.getTrain());
        Dataset<Regressor> testData= new MutableDataset<>(splitter.getTest());
        Regressor r=trainData.getExample(0).getOutput();
        System.out.println("Number of dimensions: " + r.size());
        String[] dimensionNames=r.getNames();
        System.out.println("Dimension name: " + dimensionNames[0]);
        double[] regressedValues=r.getValues();
        System.out.println("Dimension value: " + regressedValues[0]);
        Regressor.DimensionTuple tuple=r.getDimension("DIM-0").get();
        System.out.println("Tuple: ["+tuple+"]");
        LinearSGDTrainer lsgd=new LinearSGDTrainer(new SquaredLoss(), // loss function
                SGD.getLinearDecaySGD(0.01), // gradient descent algorithm
                10,                // number of training epochs
                trainData.size()/4,// logging interval
                1,                 // minibatch size
                1L                 // RNG seed
        );
        LinearSGDTrainer lrada=new LinearSGDTrainer(
                new SquaredLoss(),
                new AdaGrad(0.01),
                10,
                trainData.size()/4,
                1,
                1L
        );
        CARTRegressionTrainer cartRegressionTrainer=new CARTRegressionTrainer(6);
        XGBoostRegressionTrainer xgb=new XGBoostRegressionTrainer(50);
        Model<Regressor> lrsgdmodel=train("Linear Regression (SGD)", lsgd, trainData);
        evaluate(lrsgdmodel, testData);
        Model<Regressor> lradaModel=train("AdaGrad", lrada, trainData);
        evaluate(lradaModel, testData);
        Model<Regressor> cartModel=train("CART", cartRegressionTrainer, trainData);
        evaluate(cartModel, testData);
        Model<Regressor> xgbModel=train("XGB", xgb, trainData);
        evaluate(xgbModel, testData);


    }
    public void JDSPStuff(String path) throws IOException, WavFileException {
        Wav readTo=new Wav();
        readTo.readWav(path);
        


    }
}
