package aco_max_min;

import ecuacion.Ecuacion;

/**
 * @author ayudier
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class Ants {

    private int tourMemory[];
    private int indice;

    public Ants() {
        tourMemory = new int[1];
        tourMemory[0] = 0;
        indice = -1;
    }

    public Ants(int tamMemory) {
        this.indice = -1;
        tourMemory = new int[tamMemory];

    }

    public int newState() {
        return ++indice;
    }

    //inserta un valor(pos) al arreglo temporal
    private void inTo_Memory(int pos) {
        tourMemory[indice] = pos;
    }

    //metodo que ve las prob de los proximos estados y selecciona uno
    //modificar calcular las prob de aquellos valores a donde pueda ir
    public int nextStep(Ecuacion e, Matriz infoMedio, double pBest) {

        double sumAj = 0;
        int thisState = this.newState();
//        indice++;
        double arrayTempValue[] = new double[infoMedio.getSizeState(thisState)];

        for (int i = 0; i < infoMedio.getSizeState(thisState); i++) {
            double b = 0.0;

            b = infoMedio.getPheromone_In(thisState, i);

//            b = analizarRango(e, thisState, infoMedio, i, b);


            sumAj += b;
            arrayTempValue[i] = b;
        }

        for (int i = 0; i < arrayTempValue.length; i++) {
            arrayTempValue[i] /= sumAj;
        }
        int state = ruleta(arrayTempValue);
//        int state = desicion(arrayTempValue);
        this.inTo_Memory(state);
//        if(indice > 6){
//         int k = indice;
//        }
        return state;

    }

    //como lo pense
    public double analizarMin(Ecuacion e, int thisState, Matriz infoMedio, int i, double b) {

        double newb = b;
        boolean bcontrol = true;
        int k = 1;

        while ((e.getMinim()[thisState - k] >= 0) && bcontrol) {
            if (this.tourMemory[thisState - k] != e.getMinim()[thisState - k]) {
                bcontrol = false;
            }
            k++;
        }

        if (bcontrol && tourMemory[thisState - k] == 0) {
            if (e.getMinim()[thisState - k] == -1) {
                if (i < e.getMinim()[thisState]) {
                    newb = 0;
                }
            }
            if (e.getMinim()[thisState - k] == -2) {
                if (i > e.getMinim()[thisState]) {
                    newb = 0;
                }
            }
        }
        return newb;

    }

    public double analizarRango(Ecuacion e, int thisState, Matriz infoMedio, int i, double b) {
        double newb = b;
        if (thisState > 0) {
            newb = analizarMin(e, thisState, infoMedio, i, newb);
            newb = analizarMin(e, thisState, infoMedio, i, newb);
        }
        return newb;
    }

    public double analizarMax(Ecuacion e, int thisState, Matriz infoMedio, int i, double b) {

        double newb = b;
        boolean bcontrol = true;
        int k = 1;

        while ((e.getMaxim()[thisState - k] >= 0) && bcontrol) {
            if (this.tourMemory[thisState - k] != e.getMaxim()[thisState - k]) {
                bcontrol = false;
            }
            k++;
        }

        if (bcontrol && tourMemory[thisState - k] == infoMedio.getSizeState(thisState)) {
            if (e.getMaxim()[thisState - k] == -1) {
                if (i > e.getMaxim()[thisState]) {
                    newb = 0;
                }
            }
            if (e.getMaxim()[thisState - k] == -2) {
                if (i < e.getMaxim()[thisState]) {
                    newb = 0;
                }
            }
        }

        return newb;

    }

    private int desicion(double arrayPj[]) {
        double maximo = 0;
        int numAux = 0;
        double q = (double) Math.random();

        if (q < 0.5) {
            return (int) (Math.random() * (arrayPj.length - 1));
        } else {
            maximo = arrayPj[0];
            for (int i = 1; i < arrayPj.length; i++) {
                if (arrayPj[i] > maximo) {
                    maximo = arrayPj[i];
                    numAux = i;
                }
            }
            return numAux;
        }
    }

    private int ruleta(double arrayPj[]) {

        double minimo = 0;
        double maximo = 0;
        int numAux = 0;
        double q = (double) Math.random();

        for (int i = 0; i < arrayPj.length; i++) {
            maximo += arrayPj[i];
            if ((maximo >= q) && (q > minimo)) {
                numAux = i;
                break;
            }
            minimo = maximo;
        }
        return numAux;
    }

    public int[] returnTour() {
        return tourMemory;
    }

    public void setTourMemory(int[] tourMemory) {
        this.tourMemory = tourMemory;
    }
}
