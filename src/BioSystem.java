import java.util.ArrayList;
import java.util.Random;

public class BioSystem {

    private int L, K, s, s_max;
    private double c, alpha, timeElapsed;

    private Microhabitat[] microhabitats;
    Random rand = new Random();

    public BioSystem(int L, int S, double alpha){

        this.L = L;
        this.s = S;
        this.s_max = S;
        this.alpha = alpha;
        this.microhabitats = new Microhabitat[L];
        this.timeElapsed = 0.;

        for(int i = 0; i < L; i++){
            double c_i = Math.exp(alpha*(double)i) - 1.;
            microhabitats[i] = new Microhabitat(c_i, S);
        }
        microhabitats[0].fillWithWildType();
    }

    //this constructor is used to make BioSystems with uniform drug concentration
    public BioSystem(int L, int S, double c, char token){
        this.L = L;
        this.s = S;
        this.s_max = S;
        this.alpha = alpha;
        this.microhabitats = new Microhabitat[L];
        this.timeElapsed = 0.;

        for(int i = 0; i < L; i++){
            microhabitats[i] = new Microhabitat(c, S);
        }
        microhabitats[0].fillWithWildType();
    }

    //this constructor is used to make BioSystems with linear drug gradients
    public BioSystem(int L, int S, double minC, double maxC){

        this.L = L;
        this.s = S;
        this.s_max = S;
        this.microhabitats = new Microhabitat[L];
        this.timeElapsed = 0.;
        //calculates the linear antibiotic gradient
        double linGrad = (maxC - minC)/((double)L);

        for(int i = 0; i < L; i++){

            double c_i = linGrad*i;
            microhabitats[i] = new Microhabitat(c_i, S);
        }
        microhabitats[0].fillWithWildType();
    }


    public int getL(){
        return L;
    }
    public double getTimeElapsed(){
        return timeElapsed;
    }

    public int getCurrentPopulation(){
        int runningTotal = 0;
        for(Microhabitat m : microhabitats) {
            runningTotal += m.getN();
        }
        return runningTotal;
    }

    public int[] getSpatialDistributionArray(){
        int[] mh_pops = new int[L];
        for(int i = 0; i < L; i++){
            mh_pops[i] = microhabitats[i].getN();
        }
        return mh_pops;
    }

    public double[] getGrowthRatesArray(){
        double[] mh_gRates = new double[L];
        for(int i = 0; i < L; i++){
            mh_gRates[i] = microhabitats[i].getGrowthRate();
        }
        return mh_gRates;
    }


    public void migrate(int currentL, int bacteriumIndex){

        double direction = rand.nextDouble();

        if(direction < 0.5 && currentL < (L - 1)) {

            ArrayList<Bacteria> source = microhabitats[currentL].getPopulation();
            ArrayList<Bacteria> destination = microhabitats[currentL + 1].getPopulation();

            destination.add(source.remove(bacteriumIndex));

        }else if(direction > 0.5 && currentL > (0)){

            ArrayList<Bacteria> source = microhabitats[currentL].getPopulation();
            ArrayList<Bacteria> destination = microhabitats[currentL - 1].getPopulation();

            destination.add(source.remove(bacteriumIndex));
        }
    }

    public void die(int currentL, int bacteriumIndex){

        microhabitats[currentL].removeABacterium(bacteriumIndex);
        //if(getCurrentPopulation() == 0) populationDead = true;
    }

    public void replicate(int currentL, int bacteriumIndex){
        //a nutrient unit is consumed for every replication
        microhabitats[currentL].consumeNutrients();
        //the bacterium which is going to be replicated and its associated properties
        Bacteria parentBac = microhabitats[currentL].getBacteria(bacteriumIndex);
        int m = parentBac.getM();
        Bacteria childBac = new Bacteria(m);
        microhabitats[currentL].addABacterium(childBac);
    }



    public void performAction(){

        //selects a random bacteria from the total population
        int randomIndex = rand.nextInt(getCurrentPopulation());
        int indexCounter = 0;
        int microHabIndex = 0;
        int bacteriaIndex = 0;
        forloop:
        for(int i = 0; i < getL(); i++) {

            if((indexCounter + microhabitats[i].getN()) <= randomIndex) {

                indexCounter += microhabitats[i].getN();
                continue forloop;
            } else {
                microHabIndex = i;
                bacteriaIndex = randomIndex - indexCounter;
                break forloop;
            }
        }

        Microhabitat randMicroHab = microhabitats[microHabIndex];

        int s = randMicroHab.getS(), s_max = randMicroHab.getS_max();
        double c = randMicroHab.getC();
        Bacteria randBac = randMicroHab.getBacteria(bacteriaIndex);

        double migRate = randBac.getB();
        double deaRate = randBac.getD();
        double repliRate = randBac.replicationRate(c, s, s_max);
        double R_max = 1.2;
        double rando = rand.nextDouble()*R_max;

        if(rando < migRate) migrate(microHabIndex, bacteriaIndex);
        else if(rando >= migRate && rando < (migRate + deaRate)) die(microHabIndex, bacteriaIndex);
        else if(rando >= (migRate + deaRate) && rando < (migRate + deaRate + repliRate))
            replicate(microHabIndex, bacteriaIndex);

        timeElapsed += 1./((double) getCurrentPopulation()*R_max);
        //move this to the death() method
    }




    public static void nutrientsVsAntibioticsContourPlot(){

        int nPoints = 10, nReps = 20;
        int L = 500;
        double duration = 500.;
        String filename = "slowGrowers_nutrients_vs_antibiotic-contourPlot-FINAL";
        char token =  'c';

        ArrayList<Double> sVals = new ArrayList<Double>();
        ArrayList<Double> cVals = new ArrayList<Double>();
        ArrayList<Double> popVals = new ArrayList<Double>();

        int initS = 0, finalS = 1000;
        int sIncrement = ((finalS - initS)/nPoints);
        double initC = 0., finalC = 10.;
        double cIncrement = (finalC - initC)/(double) nPoints;


        for(double c = initC; c <= finalC; c += cIncrement) {
            cVals.add(c);

            for(int s = initS; s <= finalS; s += sIncrement) {
                sVals.add((double)s);

                double avgMaxPopulation = 0.;

                for(int r = 0; r < nReps; r++) {
                    BioSystem bs = new BioSystem(L, s, c, token);

                    while(bs.getTimeElapsed() <= duration) {
                        bs.performAction();
                    }
                    avgMaxPopulation += bs.getCurrentPopulation();
                    System.out.println(bs.getCurrentPopulation());
                    System.out.println("sVal: " + s + "\t cVal: " + c + "\t rep: " + r);
                }

                popVals.add(avgMaxPopulation/(double)nReps);
            }
        }
        System.out.println(sVals.size() +"\t"+cVals.size()+"\t"+popVals.size());
        Toolbox.writeContoursToFile(cVals, sVals, popVals, filename);
    }


    public static void nutrientsVsExponentialAntibioticsContourPlot(){

        int nPoints = 10, nReps = 20;
        int L = 500;
        double duration = 500.;
        String filename = "slowGrowers_nutrients_vs_exponential-antibiotic-contourPlot-FINAL";

        ArrayList<Double> sVals = new ArrayList<Double>();
        ArrayList<Double> alphaVals = new ArrayList<Double>();
        ArrayList<Double> popVals = new ArrayList<Double>();

        int initS = 0, finalS = 1000;
        int sIncrement = ((finalS - initS)/nPoints);
        double initAlpha = 0., finalAlpha = 0.1;
        double alphaIncrement = (finalAlpha - initAlpha)/(double)nPoints;


        for(double alpha = initAlpha; alpha <= finalAlpha; alpha += alphaIncrement) {
            alphaVals.add(alpha);

            for(int s = initS; s <= finalS; s += sIncrement) {
                sVals.add((double)s);

                double avgMaxPopulation = 0.;

                for(int r = 0; r < nReps; r++) {
                    BioSystem bs = new BioSystem(L, s, alpha);

                    while(bs.getTimeElapsed() <= duration) {
                        bs.performAction();
                    }
                    avgMaxPopulation += bs.getCurrentPopulation();
                    System.out.println(bs.getCurrentPopulation());
                    System.out.println("sVal: " + s + "\t alphaVal: " + alpha + "\t rep: " + r);
                }

                popVals.add(avgMaxPopulation/(double) nReps);
            }
        }
        System.out.println(sVals.size() +"\t"+alphaVals.size()+"\t"+popVals.size());
        Toolbox.writeContoursToFile(alphaVals, sVals, popVals, filename);
    }


    public static void linearGradient_spatialAndGRateDistributions(){

        int L = 500, nReps = 20;
        int nTimeMeasurements = 20;
        double duration = 2000., interval = duration/(double)nTimeMeasurements;
        int S = 500;
        double smallestC = 0., largestC = 12.;

        String filename = "slowGrowers-linearGradient-spatialDistribution-FINAL";
        String filename_gRate = "slowGrowers-linearGradient-gRateDistribution-FINAL";

        int[][][] allMeasurements = new int[nReps][][];
        double[][][] allGRateMeasurements = new double[nReps][][];

        for(int r = 0; r < nReps; r++){

            boolean alreadyRecorded = false;
            int[][] popsOverTime = new int[nTimeMeasurements+1][];
            double[][] gRatesOverTime = new double[nTimeMeasurements+1][];
            int timerCounter = 0;

            BioSystem bs = new BioSystem(L, S, smallestC, largestC);

            while(bs.timeElapsed <= duration){

                bs.performAction();

                if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.01) && !alreadyRecorded){

                    System.out.println("rep: "+r+"\ttime elapsed: "+String.valueOf(bs.getTimeElapsed()));
                    popsOverTime[timerCounter] = bs.getSpatialDistributionArray();
                    gRatesOverTime[timerCounter] = bs.getGrowthRatesArray();

                    alreadyRecorded = true;
                    timerCounter++;
                }
                if(bs.getTimeElapsed()%interval >= 0.1) alreadyRecorded = false;
            }

            allMeasurements[r] = popsOverTime;
            allGRateMeasurements[r] = gRatesOverTime;
        }

        double[][] averagedPopDistributions = Toolbox.averagedResults(allMeasurements);
        double[][] averagedGRateDistributions = Toolbox.averagedResults(allGRateMeasurements);
        Toolbox.printAveragedResultsToFile(filename, averagedPopDistributions);
        Toolbox.printAveragedResultsToFile(filename_gRate, averagedGRateDistributions);
    }


    public static void exponentialGradient_spatialAndGRateDistributions(double input_alpha){

        int L = 500, nReps = 20;
        int nTimeMeasurements = 20;
        double duration = 2000., interval = duration/(double)nTimeMeasurements;
        double alpha = input_alpha;
        int S = 500;

        String filename = "slowGrowers-alpha="+String.valueOf(alpha)+"-spatialDistribution-FINAL";
        String filename_gRate = "slowGrowers-alpha="+String.valueOf(alpha)+"-gRateDistribution-FINAL";

        int[][][] allMeasurements = new int[nReps][][];
        double[][][] allGRateMeasurements = new double[nReps][][];

        for(int r = 0; r < nReps; r++){

            boolean alreadyRecorded = false;
            int[][] popsOverTime = new int[nTimeMeasurements+1][];
            double[][] gRatesOverTime = new double[nTimeMeasurements+1][];
            int timerCounter = 0;

            BioSystem bs = new BioSystem(L, S, alpha);

            while(bs.timeElapsed <= duration){

                bs.performAction();

                if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.01) && !alreadyRecorded){

                    System.out.println("rep: "+r+"\ttime elapsed: "+String.valueOf(bs.getTimeElapsed()));
                    popsOverTime[timerCounter] = bs.getSpatialDistributionArray();
                    gRatesOverTime[timerCounter] = bs.getGrowthRatesArray();

                    alreadyRecorded = true;
                    timerCounter++;
                }
                if(bs.getTimeElapsed()%interval >= 0.1) alreadyRecorded = false;
            }

            allMeasurements[r] = popsOverTime;
            allGRateMeasurements[r] = gRatesOverTime;
        }

        double[][] averagedPopDistributions = Toolbox.averagedResults(allMeasurements);
        double[][] averagedGRateDistributions = Toolbox.averagedResults(allGRateMeasurements);
        Toolbox.printAveragedResultsToFile(filename, averagedPopDistributions);
        Toolbox.printAveragedResultsToFile(filename_gRate, averagedGRateDistributions);
    }
}

























