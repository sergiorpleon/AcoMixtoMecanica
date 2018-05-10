/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ecuacion;

import java.io.*;
import mecanica.*;

/**
 *
 * @author usuario
 */
public class Ecuacion {
    
    private int gl;
    private int cant_discreta;
    private int cant_continua;
    private double[] min;       //para variante de acotar
    private double[] max;       //para variante de acotar
    private int[] minim = new int[0];       //para variante de seleccion
    private int[] maxim = new int[0];       //para variante de seleccion
    private String[] nombre;
    private String[] significado;
    private String[] type;
    private String[][] dominio;
    private int indice;
    
    private Intercambiadores itc;

    public Ecuacion(int gl, int cant_discreta, int cant_continua, String[] nombre, String[] significado, String[] type, String[][] dominio, double[] min, double[] max) {
        this.gl = gl;
        this.cant_discreta = cant_discreta;
        this.cant_continua = cant_continua;
        this.nombre = nombre;
        this.significado = significado;
        this.type = type;
        this.dominio = dominio;

        this.min = min;
        this.max = max;
        this.minim = new int[cant_discreta + cant_continua * (gl + 1)];
        this.maxim = new int[cant_discreta + cant_continua * (gl + 1)];

        indice = cant_continua + cant_discreta;
        
    }

    public Ecuacion(int cant_variable, int gl) {
        this.gl = gl;
        this.nombre = new String[cant_variable];
        this.significado = new String[cant_variable];
        this.type = new String[cant_variable];
        this.dominio = new String[cant_variable][];
        this.min = new double[0];
        this.max = new double[0];
        indice = 0;
        
    }

    /*
     *crea los arreglos min max para el metodo de actoar
     */
    private void addMinMax(double newmin, double newmax) {
        double[] tmpmin = new double[cant_continua];
        double[] tmpmax = new double[cant_continua];
        for (int i = 0; i < cant_continua - 1; i++) {
            tmpmin[i] = min[i];
            tmpmax[i] = max[i];
        }
        tmpmin[cant_continua - 1] = newmin;
        tmpmax[cant_continua - 1] = newmax;

        min = tmpmin;
        max = tmpmax;
    }

    /*
     *crean los arreglos minim y maxim para el metodo de seleccion
     */
    private void addMinimMaxim() {
        int[] tmpminim = new int[cant_discreta + cant_continua * (gl + 1)];
        int[] tmpmaxim = new int[cant_discreta + cant_continua * (gl + 1)];
        for (int i = 0; i < this.minim.length; i++) {
            tmpminim[i] = minim[i];
            tmpmaxim[i] = maxim[i];
        }
        tmpminim[this.minim.length] = -1;
        tmpmaxim[this.minim.length] = -1;

        minim = tmpminim;
        maxim = tmpmaxim;

    }

    private void addMinimMaxim(double newmin, double newmax) {
        int[] tmpminim = new int[cant_discreta + cant_continua * (gl + 1)];
        int[] tmpmaxim = new int[cant_discreta + cant_continua * (gl + 1)];
        for (int i = 0; i < this.minim.length; i++) {
            tmpminim[i] = minim[i];
            tmpmaxim[i] = maxim[i];
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int sigmin = 1;
        int sigmax = 1;

        if (newmin > 0) {
            tmpminim[minim.length] = -1;
        } else {
            tmpminim[minim.length] = -2;
            sigmin = -1;
        }
        if (newmax > 0) {
            tmpmaxim[minim.length] = -1;
        } else {
            tmpmaxim[minim.length] = -2;
            sigmax = -1;
        }

        int nummin = (int) (newmin * Math.pow(10, gl));
        int nummax = (int) (newmax * Math.pow(10, gl));

        for (int k = 2; k > 0; k--) {
            tmpminim[minim.length + k] = (sigmin * nummin) % 10;
            nummin = (nummin - ((sigmin * nummin) % 10)) / 10;

            tmpmaxim[minim.length + k] = (sigmax * nummax) % 10;
            nummax = (nummax - ((sigmin * nummax) % 10)) / 10;
        }

        minim = tmpminim;
        maxim = tmpmaxim;
    }

    public void add_variable(String nombre, String significado, String type, double min, double max) {
        this.nombre[indice] = nombre;
        this.significado[indice] = significado;
        this.type[indice] = type;
        int intmin = 0;
        int intmax = 0;

        if (min < 0) {
            intmin = (int) (Math.floor(min) + 1);
        } else {
            intmin = (int) Math.floor(min);
        }

        if (max < 0) {
            intmax = (int) (Math.floor(max) + 1);
        } else {
            intmax = (int) Math.floor(max);
        }


        this.dominio[indice] = new String[intmax - intmin + 1];
        for (int i = 0; i < this.dominio[indice].length; i++) {

            this.dominio[indice][i] = intmin + i + "";

        }


        /*caso -0*/
        if (this.type[indice].compareToIgnoreCase("continua") == 0) {
            if (dominio[indice][0].compareToIgnoreCase("0") == 0) {
                /*caso -0 a x*/
                if (min < 0) {
                    //insertar -0 al inicio
                    String[] tmpdominio = new String[dominio[indice].length + 1];
                    for (int i = 1; i < tmpdominio.length; i++) {
                        tmpdominio[i] = dominio[indice][i - 1];
                    }
                    tmpdominio[0] = "-0";

                    this.dominio[indice] = new String[tmpdominio.length];
                    for (int i = 0; i < tmpdominio.length; i++) {

                        this.dominio[indice][i] = tmpdominio[i];

                    }

                }
            } else if (dominio[indice][dominio[indice].length - 1].compareToIgnoreCase("0") == 0) {
                /*caso -x a -0*/
                if (max < 0) {
                    //cambiar 0 por -0
                    String[] tmpdominio = new String[dominio[indice].length];
                    for (int i = 0; i < dominio[indice].length; i++) {
                        tmpdominio[i] = dominio[indice][i];
                    }
                    tmpdominio[dominio[indice].length - 1] = "-0";

                    this.dominio[indice] = new String[tmpdominio.length];
                    for (int i = 0; i < this.dominio[indice].length; i++) {

                        this.dominio[indice][i] = tmpdominio[i];

                    }
                    /*caso -x a 0*/
                } else {
                    //adicionar -0 en penultima
                    String[] tmpdominio = new String[dominio[indice].length + 1];
                    for (int i = 0; i < dominio.length; i++) {
                        tmpdominio[i] = dominio[indice][i];
                    }
                    tmpdominio[dominio.length - 1] = "-0";
                    tmpdominio[dominio.length] = "0";

                    this.dominio[indice] = new String[tmpdominio.length];
                    for (int i = 0; i < this.dominio[indice].length; i++) {

                        this.dominio[indice][i] = tmpdominio[i];

                    }
                }
            } else /*caso -x x*/ if (min * max < 0) {
                //adicionar
                String[] tmpdominio = new String[dominio[indice].length + 1];
                int i = 0;
                while (dominio[indice][i].compareToIgnoreCase("0") != 0) {
                    tmpdominio[i] = dominio[indice][i];
                    i++;
                }
                tmpdominio[i++] = "-0";

                while (i < tmpdominio.length) {
                    tmpdominio[i] = dominio[indice][i - 1];
                    i++;
                }

                this.dominio[indice] = new String[tmpdominio.length];
                for (int j = 0; j < this.dominio[indice].length; j++) {

                    this.dominio[indice][j] = tmpdominio[j];

                }
            }

        }
        /*caso -0*/

        if (this.type[indice].compareToIgnoreCase("discreta") == 0) {
            this.cant_discreta++;

            addMinimMaxim();

        } else {
            this.cant_continua++;
            addMinMax(min, max);

            addMinimMaxim(min, max);
        }
        indice++;


    }

    public void add_variable(String nombre, String significado, String type, String[] dominio) {
        this.nombre[indice] = nombre;
        this.significado[indice] = significado;
        this.type[indice] = type;
        this.dominio[indice] = dominio;
        this.cant_discreta++;
        addMinimMaxim();
        indice++;
    }

    public double eval_funcion_mecanica(double[] variable) {
      //public Intercambiadores             ( fluido_masico_coraza,  fluido_masico_tubo,  calor_especifico_caliente,  calor_especifico_frio, double temperatura_entrada_caliente,  temperatura_salida_caliente,  temperatura_entrada_frio,  temperatura_salda_frio,  resistencia_ensuciamiento_tubo,  resistencia_ensuciamiento_coraza,  conductividad_termica_tubo,  conductividad_termica_coraza,  densidad_tubo,  densidad_coraza,   viscosidad_dinamica_tubo,  viscosidad_dinamica_coraza,  numero_annos_vida,  tasa_descuento_anual,  costo_energia_electrica,  horas_trabajo_anno) {
        return itc.evaluar(variable);
    }
    
    public String[] imprimir_mecanica(double[] variable) {
      //public Intercambiadores             ( fluido_masico_coraza,  fluido_masico_tubo,  calor_especifico_caliente,  calor_especifico_frio, double temperatura_entrada_caliente,  temperatura_salida_caliente,  temperatura_entrada_frio,  temperatura_salda_frio,  resistencia_ensuciamiento_tubo,  resistencia_ensuciamiento_coraza,  conductividad_termica_tubo,  conductividad_termica_coraza,  densidad_tubo,  densidad_coraza,    viscosidad_dinamica_tubo,  viscosidad_dinamica_coraza,  numero_annos_vida,  tasa_descuento_anual,  costo_energia_electrica,  horas_trabajo_anno) {
        return itc.imprimir(variable);
    }
    
//    public double eval_funcion(double[] variable, int alpha, int beta) {
//        double evaluacion = 0;
//        int[] var_discreta = new int[cant_discreta];
//        double[] var_continua = new double[cant_continua];
//        int cont_discreta = 0;
//        int cont_continua = 0;
//
//        for (int i = 0; i < variable.length; i++) {
//            if (this.type[i].compareToIgnoreCase("discreta") == 0) {
//                var_discreta[cont_discreta] = (int) variable[i];
//                cont_discreta++;
//            } else {
//                var_continua[cont_continua] = variable[i];
//                cont_continua++;
//            }
//
//        }
//
////        evaluacion = funcion.evaluar(var_discreta, var_continua);
//
//        evaluacion = funcion.evaluarDiscreta(alpha, var_discreta);
//
//        evaluacion += bf.eval(beta, var_continua);
//
//        return evaluacion;
//
//    }

//    public int signo()
    public int getSize() {
        return (cant_continua * (gl + 1) + cant_discreta);
    }

    public double[][] construir_matriz(double phermoneMaximum) {
        int fila = 0;
        int num_fila = cant_continua * (gl + 1) + cant_discreta;
        double matriz[][];
        matriz = new double[num_fila][];


        for (int i = 0; i < cant_continua + cant_discreta; i++) {
            matriz[fila] = new double[dominio[i].length];

            for (int j = 0; j < dominio[i].length; j++) {
                matriz[fila][j] = phermoneMaximum;
            }


            if (type[i].compareToIgnoreCase("continua") == 0) {

                for (int k = 0; k < gl; k++) {
                    matriz[fila + k + 1] = new double[10];

                    for (int j = 0; j < 10; j++) {
                        matriz[fila + k + 1][j] = phermoneMaximum;
                    }
                }
                fila += gl;
            }
            fila++;
        }

//        for (int i = 0; i < num_fila; i++) {
//            for (int j = 0; j < matriz[i].length; j++) {
//                System.out.print(matriz[i][j] + "\t");
//            }
//            System.out.println();
//        }


        return matriz;
    }

    public double[] construir_variable(int[] tour) {
        int fila = 0;
        double[] vector = new double[cant_continua + cant_discreta];
        double num = 0;
        int signo;
        for (int i = 0; i < vector.length; i++) {
            vector[i] = Double.parseDouble(dominio[i][tour[fila]]);

            if (type[i].compareToIgnoreCase("continua") == 0) {
                signo = 1;
                //conformo el numero
                num = 0;
                //tomo el valor que refleja tour[fila en el dominio]
                num += Double.parseDouble(this.dominio[i][tour[fila]]);
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if (Double.parseDouble(this.dominio[i][tour[fila]]) < 0 || this.dominio[i][tour[fila]].compareToIgnoreCase("-0") == 0) {
                    signo = -1;
                }
                for (int k = 0; k < gl; k++) {
                    num += signo * tour[fila + k + 1] / (Math.pow(10, k + 1));
                }

                num = Math.round(num * Math.pow(10, gl)) / Math.pow(10, gl);

                vector[i] = num;
                fila += gl;
            }
            fila++;
        }
        return vector;
    }

    public int[] getMaxim() {
        return maxim;
    }

    public int[] getMinim() {
        return minim;
    }

    public double[] getMax() {
        return max;
    }

    public double[] getMin() {
        return min;
    }

    public int getCant_contontinua() {
        return cant_continua;
    }

    public void setCant_contontinua(int cant_contontinua) {
        this.cant_continua = cant_contontinua;
    }

    public int getCant_discreta() {
        return cant_discreta;
    }

    public void setCant_discreta(int cant_discreta) {
        this.cant_discreta = cant_discreta;
    }

    public String[][] getDominio() {
        return dominio;
    }

    public void setDominio(String[][] dominio) {
        this.dominio = dominio;
    }

    public int getGl() {
        return gl;
    }

    public void setGl(int gl) {
        this.gl = gl;
    }

    public String[] getSignificado() {
        return significado;
    }

    public void setSignificado(String[] significado) {
        this.significado = significado;
    }

    public String[] getType() {
        return type;
    }

    public void setType(String[] type) {
        this.type = type;
    }

    public String[] getVariable() {
        return nombre;
    }

    public void setVariable(String[] variable) {
        this.nombre = variable;
    }

    public Ecuacion(String filePath, Intercambiadores itc) {
        this.itc = itc;
        try {
            int textgl;
            int textcant_variable;
            double textmin;
            double textmax;
            String textnombre;
            String textdefinicion;
            String texttype;
            String[] textdominio;


            File file = new File(filePath);
            FileReader f = new FileReader(file);
            String s = new String();
            BufferedReader a = new BufferedReader(f);
            s = a.readLine();
            s = s.split(": ")[1];
            textcant_variable = Integer.parseInt(s);
            s = a.readLine();
            s = s.split(": ")[1];
            textgl = Integer.parseInt(s);

            /*poner en una funcion
             */
            this.gl = textgl;
            this.nombre = new String[textcant_variable];
            this.significado = new String[textcant_variable];
            this.type = new String[textcant_variable];
            this.dominio = new String[textcant_variable][];
            this.min = new double[0];
            this.max = new double[0];
            indice = 0;
            
            /*
             */

            String[] cadArray1;
            String[] cadArray2;
            String[] cadArray3;
            String[] cadArray4;
            String[] cadArray5;
            String cad1;

            for (int i = 0; i < textcant_variable; i++) {
                s = a.readLine();

                cadArray1 = s.split("[\"]+"); // nombre definicion resto
                cadArray2 = cadArray1[0].split("[ \t]+"); //nombre
                cadArray3 = cadArray1[2].split("[{}]+"); //tipo-posmaxmin dominio
                cadArray4 = cadArray3[0].split("[ \t]+"); //tipo


                textnombre = cadArray2[0];
                textdefinicion = cadArray1[1];
                texttype = cadArray4[1];

                if (cadArray3.length == 1) {
                    textmin = Double.parseDouble(cadArray4[2]);
                    textmax = Double.parseDouble(cadArray4[3]);

                    this.add_variable(textnombre, textdefinicion, texttype, textmin, textmax);

                } else {
                    cad1 = cadArray3[1];
                    cadArray5 = cad1.split("[, \t]+");
                    textdominio = cadArray5;

                    this.add_variable(textnombre, textdefinicion, texttype, textdominio);
                }
            }
            a.close();

        } catch (Exception e) {
            System.out.println("Ooops.  Error: " + e.getMessage());
        }
    }
//    con descripcion al final
//    public void readFile(String filePath) {
//        try {
//            int gl;
//            int cant_variable;
//            int min;
//            int max;
//            String nombre;
//            String definicion;
//            String type;
//            String[] dominio;
//
//
//            File file = new File(filePath);
//            FileReader f = new FileReader(file);
//            String s = new String();
//            BufferedReader a = new BufferedReader(f);
//            s = a.readLine();
//            cant_variable = Integer.parseInt(s);
//            s = a.readLine();
//            gl = Integer.parseInt(s);
//
//            ecuacion = new Ecuacion(cant_variable, gl);
//
//            String[] cadArray1;
//            String[] cadArray2;
//            String[] cadArray3;
//            String[] cadArray4;
//            String[] cadArrayaux;
//            String cad1;
//            String cad2;
//
//            for (int i = 0; i < cant_variable; i++) {
//                s = a.readLine();
//
//                cadArray1 = s.split("[\"]+");
//                cadArray2 = cadArray1[0].split("[{}]+");
//
//
//                cadArray3 = cadArray2[0].split("[ \t]+");
//
//
//                nombre = cadArray3[0];
//                definicion = cadArray1[1];
//                type = cadArray3[1];
//
//                if (cadArray2.length == 1) {
//                    min = Integer.parseInt(cadArray3[2]);
//                    max = Integer.parseInt(cadArray3[3]);
//
//                    ecuacion.add_variable(nombre, definicion, type, min, max);
//                } else {
//                    cad2 = cadArray2[1];
//                    cadArray4 = cad2.split("[; \t]+");
//                    dominio = cadArray4;
//
//                    ecuacion.add_variable(nombre, definicion, type, dominio);
//                }
//            }
//            a.close();
//
//        } catch (Exception e) {
//            System.out.println("Ooops.  Error: " + e.getMessage());
//        }
//    }
}
