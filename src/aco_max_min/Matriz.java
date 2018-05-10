/**
 * 
 */
package aco_max_min;

import ecuacion.Ecuacion;

/**
 * @author puris Jun 6, 2006 
 *
 */
public class Matriz {

    private double matrixPheromone[][];
    private double pheromoneInicial;
    private double pheromoneMaximum;
    private double pheromoneMinimum;
    private double constEvaporation;
    private int sizeProblem;

    /**
     *en esta clase se realiza la mayoria de las transformaciones (adaptaciones del problema)
     * no se cuenta con la matriz de peso en su lugar se utliza la evaluacion en una funcion
     * se debe crear la matriz de fermona
     */
    public Matriz(Ecuacion e, double constEvaporation) {
        matrixPheromone = e.construir_matriz(pheromoneMaximum);
        sizeProblem = e.getSize();
        this.pheromoneMaximum = constEvaporation * 1 / 1000 + constEvaporation;
        this.pheromoneMinimum = this.pheromoneMaximum / (2 * this.sizeProblem);
        this.pheromoneInicial = this.pheromoneMaximum;
        this.constEvaporation = constEvaporation;
        this.inicializeMatrixPheromone();

    }

    public void printMaxMinPheromone() {
        System.out.println("Maximum : " + this.pheromoneMaximum);
        System.out.println("Minimum : " + this.pheromoneMinimum);
    }

    public void printMatrixPheromone() {
        for (int i = 0; i < this.matrixPheromone.length; i++) {
            for (int j = 0; j < this.matrixPheromone[i].length; j++) {
//                System.out.print((Math.floor(this.matrixPheromone[i][j]*100000))/100000 + " ");
                System.out.print(this.matrixPheromone[i][j] + " ");
            }
            System.out.println();
        }
    }

    // este metodo calcula las cuotas superior e inferior de la feromona.
    /**
     * @param cost :el valor de la mejor solucion en la iteracion.
     * @param pBest : probabilidad de eleccion de cada candidato.
     * @param iteration: la iteracion actual.
     */
    public void calculate_PheromoneMaxMin(double cost, double pBest) {
        this.pheromoneMaximum = (1 / (1 - this.constEvaporation)) * (1 / cost);
        this.pheromoneMinimum = this.pheromoneMaximum / (2.0 * this.sizeProblem);
    }

    // este metodo es el encargado de ajustar los valores de feromana para el Max Min
    public void regulate_PheromoneMaxMin() {
        for (int i = 0; i < this.matrixPheromone.length; i++) {
            for (int j = 0; j < this.matrixPheromone[i].length; j++) {
                if (this.matrixPheromone[i][j] > this.pheromoneMaximum) {
                    this.matrixPheromone[i][j] = this.pheromoneMaximum;
                } else if (this.matrixPheromone[i][j] < this.pheromoneMinimum) {
                    this.matrixPheromone[i][j] = this.pheromoneMinimum;
                }
            }
        }
    }

    /**
     * @param tour se realiza el update de la ferromona global y primero se hace una evaporacion
     */
    // funciones de update de ferromona para el algoritmo AS
    public void updatePheromone_Global(double coste, int tour[]) {
        for (int i = 0; i < this.matrixPheromone.length; i++) {
            for (int j = 0; j < this.matrixPheromone[i].length; j++) {
                this.matrixPheromone[i][j] = (1 - this.constEvaporation) * this.matrixPheromone[i][j];
            }
        }
        //ver con detenimiento, actualiza pheromona en el camino
        for (int i = 0; i < tour.length; i++) {
            this.matrixPheromone[i][tour[i]] += this.constEvaporation * (1 / coste);
        }

    }

    /**simbolizan el paso que dio la hormiga, se hace el update paso a paso.
     */
    public void updatePheromone_Step_by_Step(int i, int j) {
        this.matrixPheromone[i][j] = (1 - this.constEvaporation) * this.matrixPheromone[i][j] + this.constEvaporation * this.pheromoneInicial;

    }

    public double getPheromone_In(int i, int j) {
        return this.matrixPheromone[i][j];
    }

//
    public int getSizeState(int i) {
        return matrixPheromone[i].length;
    }

    private void inicializeMatrixPheromone() {
        for (int i = 0; i < this.matrixPheromone.length; i++) {
            for (int j = 0; j < this.matrixPheromone[i].length; j++) {
                this.matrixPheromone[i][j] = this.pheromoneMaximum;
            }
        }
    }
    //cambiaria completamente
}
