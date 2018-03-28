import java.util.ArrayList;

public class Microhabitat {

    private int s, s_max, K = 100;
    private double c;

    private final double K_prime = 33.;

    private ArrayList<Bacteria> population;

    public Microhabitat(double c, int S){
        this.c = c;
        this.s = S;
        this.s_max = S;
        this.population = new ArrayList<Bacteria>();
    }

    public double getC(){return c;}
    public int getS(){return s;}
    public int getS_max(){return s_max;}

    public int getN(){
        return population.size();
    }

    public ArrayList<Bacteria> getPopulation(){
        return population;
    }
    public Bacteria getBacteria(int i){
        return population.get(i);
    }

    public void fillWithWildType(){
        int initGenotype = 1;

        for(int i = 0; i < K; i++){
            population.add(new Bacteria(initGenotype));
        }
    }

    public double getGrowthRate(){
        double mu = s/(K_prime+s);
        double mu_max = s_max/(K_prime+s_max);
        double beta = 1. + 9.*(mu/mu_max);
        double phi_c = 1. - (c/beta)*(c/beta);

        return (phi_c > 0.) ? phi_c*(s/(K_prime+s)) : 0.;
    }

    public void consumeNutrients(){
        if(s > 0) s--;
    }

    public void removeABacterium(int i){
        population.remove(i);
    }

    public void addABacterium(Bacteria newBact){
        population.add(newBact);
    }
}
