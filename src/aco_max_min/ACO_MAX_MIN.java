/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package aco_max_min;

import java.util.Date;
import ecuacion.Ecuacion;
import java.awt.image.BufferedImage;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;
import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.data.*;
import org.jfree.data.xy.*;



/**
 *
 * @author usuario
 */
public class ACO_MAX_MIN {

    private Ecuacion ecuacion;
    private Ants arrayAnts[];
    private Matriz pheromoneMatriz;
    private int[] bestGlobal_Solution;
    private int cantCicle;
    private int cantAnts;
    private double pBest;
    private int iteration;

    public ACO_MAX_MIN(double pBest, double constEvaporacion, int cantAnts, int cantCicle, Ecuacion e1) {
        this.ecuacion = e1;
        this.pheromoneMatriz = new Matriz(ecuacion, constEvaporacion);
        this.cantAnts = cantAnts;
        this.cantCicle = cantCicle;
        this.pBest = pBest;
    }

    public ACO_MAX_MIN(double constEvaporacion, int cantAnts, int cantCicle, Ecuacion e1) {
        this.ecuacion = e1;
        this.pheromoneMatriz = new Matriz(ecuacion, constEvaporacion);
        this.cantAnts = cantAnts;
        this.cantCicle = cantCicle;

    }

    //seria evaluar la mejor solucion y copararla con la mejor global
    //calcular coste globar, debe ser el metodo eval function, en la clase medio
    private boolean CompareGlobalSolution(int bestCycle_Solution[]) {
        double sumSolutioAux = 0;
        double sumGlobalSolutio = 0;
        if (this.bestGlobal_Solution != null) {
            sumSolutioAux = calculateCoste_Global(bestCycle_Solution);        //seria eval function con solution mejor
            sumGlobalSolutio = calculateCoste_Global(this.bestGlobal_Solution);//seria eval funtion con solution mejor global
            if (sumGlobalSolutio > sumSolutioAux) {
                return true;
            }
            return false;
        } else {
            return true;
        }
    }

    private void PrintGlobalSolution() {
        System.out.print("Costo : ");
        System.out.println(calculateCoste_Global(this.bestGlobal_Solution));
        System.out.println("Iteracion : " + this.iteration);
    }

    //busqueda del optimo
    private int[] localSerch2_opt(int[] soluction) {
//
        int[] arraytemp = new int[soluction.length];
        double soluction_result;

        for (int i = 0; i < soluction.length; i++) {
            arraytemp[i] = soluction[i];
        }
        soluction_result = calculateCoste_Global(soluction);

        for (int i = 0; i < this.arrayAnts.length; i++) {
            for (int j = 0; j < soluction.length; j++) {
                arraytemp[j] = this.arrayAnts[i].returnTour()[j];
                /*1*/
//                   arraytemp = this.informationMedio.acotar(arraytemp); //acotarrrrrrrrrrrrrr
                   /*1*/
                if (esAcotado(arraytemp)) {
                    if (soluction_result > calculateCoste_Global(arraytemp)) {
                        soluction[j] = arraytemp[j];
                        soluction_result = calculateCoste_Global(soluction);
                    } else {
                        arraytemp[j] = soluction[j];
                    }
                } else {
                    arraytemp[j] = soluction[j];
                }
            }
        }
        return arraytemp;
    }

    public void executeAsTSP(JLabel jl) {
        XYSeries series = null;
        XYSeries series2 = null;
        XYDataset datos;
        JFreeChart linea = null;
        series= new XYSeries("Valor de la mejor hormiga de la iteración");
        series2= new XYSeries("Valor de la mejor hormiga global");
        BufferedImage graficoLinea = null;


        
        
        this.arrayAnts = new Ants[this.cantAnts];
        int cantState = getSizeProblem();
        this.pheromoneMatriz.printMaxMinPheromone();

        for (int i = 0; i < cantCicle; i++) { // cantidad de ciclos (iteraciones)
            for (int j = 0; j < this.cantAnts; j++) { // inicializacion del arreglo de hormigas y del estado inicial
                arrayAnts[j] = new Ants(cantState);
                arrayAnts[j].nextStep(ecuacion, this.pheromoneMatriz, this.pBest);
            }

            int k = 1;
            //recorro # de fila
            while (k++ < cantState) {
                for (int j = 0; j < this.cantAnts; j++) {// movimiento de las hormigas
                    arrayAnts[j].nextStep(ecuacion, this.pheromoneMatriz, this.pBest);
                }
            }

            // find a local solution

            /*1*/
            for (int j = 0; j < cantAnts; j++) {
                arrayAnts[j].setTourMemory(this.acotar(arrayAnts[j].returnTour()));   //acotarrrrrr
            }
            /*1*/

            int[] localSolution = this.arrayAnts[0].returnTour();
            for (int j = 1; j < this.cantAnts; j++) {
                if (calculateCoste_Global(localSolution) > calculateCoste_Global(this.arrayAnts[j].returnTour())) {
                    localSolution = this.arrayAnts[j].returnTour();
                }
            }

            localSolution = localSerch2_opt(localSolution);
//            localSolution = this.localSerch2_opt(localSolution);
            if (this.CompareGlobalSolution(localSolution)) {
                this.bestGlobal_Solution = localSolution;
                this.iteration = i;
            }

            double coste = calculateCoste_Global(localSolution);
            this.pheromoneMatriz.updatePheromone_Global(coste, localSolution);
            this.pheromoneMatriz.calculate_PheromoneMaxMin(coste, this.pBest);
            this.pheromoneMatriz.regulate_PheromoneMaxMin();
            
            
            series.add(i,coste);
            series2.add(i, calculateCoste_Global(this.bestGlobal_Solution));
        //jfreee chart
            
        XYSeriesCollection allserie = new XYSeriesCollection();
        allserie.addSeries(series);
        allserie.addSeries(series2);
        
        datos = allserie;
        
        
        linea = ChartFactory.createXYLineChart("Visualización del proceso de optimización","Iteraciones","Costo",datos,PlotOrientation.VERTICAL,true,true,true);
        
        graficoLinea=linea.createBufferedImage(400, 400);
       
        jl.setSize(400, 400);
        jl.setIcon(new ImageIcon(graficoLinea));
        jl.paint(jl.getGraphics());
        



        
        
        //
            
            
            
            
        }

        this.arrayAnts = new Ants[0];
        
        

    }

    public int getSizeProblem() {
        return ecuacion.getSize();
    }

    public double calculateCoste_Global(int tour[]) {
        double cost = 0;
        //crear vector sol
        double[] variable = ecuacion.construir_variable(tour);
        cost = ecuacion.eval_funcion_mecanica(variable);
//        cost = ecuacion.eval_funcion(variable, alpha, beta);
        return cost;
    }

    public String[] imprimir_Global() {

        //crear vector sol
        double[] variable = ecuacion.construir_variable(bestGlobal_Solution);
        return ecuacion.imprimir_mecanica(variable);
//        cost = ecuacion.eval_funcion(variable, alpha, beta);

    }

    public boolean esAcotado(int tour[]) {
        boolean b = true;
        int[] newtour = new int[tour.length];
        for (int i = 0; i < tour.length; i++) {
            newtour[i] = tour[i];
        }
        int fila = 0;
        double num = 0;
        int count = -1; //cuenta la cant de variables continuas


        //conformo el numero para compararlo com max, min
        for (int i = 0; i < ecuacion.getVariable().length; i++) {
            if (ecuacion.getType()[i].compareToIgnoreCase("continua") == 0) {
                count++;
                //conformo el numero
                num = 0;
                //tomo el valor que refleja tour[fila en el dominio]
                num += Double.parseDouble(ecuacion.getDominio()[i][tour[fila]]);
                for (int k = 0; k < ecuacion.getGl(); k++) {
                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    if (Double.parseDouble(ecuacion.getDominio()[i][tour[fila]]) < 0 || ecuacion.getDominio()[i][tour[fila]].compareToIgnoreCase("-0") == 0) {
                        num += -1 * tour[fila + k + 1] / (Math.pow(10, k + 1));
                    } else {
                        num += tour[fila + k + 1] / (Math.pow(10, k + 1));
                    }
                }

                num = Math.round(num * Math.pow(10, ecuacion.getGl())) / Math.pow(10, ecuacion.getGl());

                //verifico si no esta en el rango
                if (num < ecuacion.getMin()[count] || num > ecuacion.getMax()[count]) {
                    b = b && false;
                }
                fila += ecuacion.getGl();
            }
            fila++;
        }
        return b;
    }

    public int[] acotar(int tour[]) {
        int[] newtour = new int[tour.length];
        for (int i = 0; i < tour.length; i++) {
            newtour[i] = tour[i];
        }
        int fila = 0;
        double num = 0;
        int count = -1; //cuenta la cant de variables continuas
        int tmpnum = 0;
        int signo = 1;

        //conformo el numero para compararlo com max, min
        for (int i = 0; i < ecuacion.getVariable().length; i++) {
            if (ecuacion.getType()[i].compareToIgnoreCase("continua") == 0) {
                count++;
                //conformo el numero
                num = 0;
                //tomo el valor que refleja tour[fila en el dominio]
                num += Double.parseDouble(ecuacion.getDominio()[i][tour[fila]]);
                for (int k = 0; k < ecuacion.getGl(); k++) {
                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    if (Double.parseDouble(ecuacion.getDominio()[i][tour[fila]]) < 0 || ecuacion.getDominio()[i][tour[fila]].compareToIgnoreCase("-0") == 0) {
                        num += -1 * tour[fila + k + 1] / (Math.pow(10, k + 1));
                    } else {
                        num += tour[fila + k + 1] / (Math.pow(10, k + 1));
                    }
                }

                num = Math.round(num * Math.pow(10, ecuacion.getGl())) / Math.pow(10, ecuacion.getGl());

                //verifico si no esta en el rango
                if (num < ecuacion.getMin()[count] || num > ecuacion.getMax()[count]) {
                    signo = 1;
                    if (num < ecuacion.getMin()[count]) {
                        if (ecuacion.getMin()[count] < 0) {
                            signo = -1;
                        }
                        tmpnum = (int) (ecuacion.getMin()[count] * Math.pow(10, ecuacion.getGl()));
                    } else {
                        if (ecuacion.getMax()[count] < 0) {
                            signo = -1;
                        }
                        tmpnum = (int) (ecuacion.getMax()[count] * Math.pow(10, ecuacion.getGl()));
                    }

                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//                    if (tmpnum < 0) {
//                        signo = -1;
//                    }

                    //inserto max o min
                    tmpnum = signo * tmpnum;
                    for (int k = ecuacion.getGl(); k > 0; k--) {
                        newtour[fila + k] = tmpnum % 10;
                        tmpnum = (tmpnum - (tmpnum % 10)) / 10;
                    }


                    //el valor se indexa segun el dominio
//                    if (num < ecuacion.getMin()[count]) {
//                        newtour[fila] = 0;
//                    } else {
//                        newtour[fila] = ecuacion.getDominio()[i].length - 1;
//                    }

                }
                fila += ecuacion.getGl();
            }
            fila++;
        }
        return newtour;
    }
}
