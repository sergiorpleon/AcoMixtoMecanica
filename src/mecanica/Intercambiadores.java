/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package mecanica;

/**
 *
 * @author sperez
 */
public class Intercambiadores {
    /*
     * constantes
     */

    private static final double a1 = 8000;
    private static final double a2 = 259.2;
    private static final double a3 = 0.91;
    private static final double eficiencia_bombeo = 0.7;
    private static final double contante_p = 2.5;

    /*
     * datos
     */
    private double m_s_h;       //fluido caliente
    private double m_t_c;       //fluido frio
    private double cp_s_h;
    private double cp_t_c;
    private double Ti_s_h;
    private double To_s_h;
    private double Ti_t_c;
    private double To_t_c;
    private double Rensuc_t;
    private double Rensuc_s;
    private double LANDA_t;
    private double LANDA_s;
    private double RO_t;
    private double RO_s;
    private double NIU_t;
    private double NIU_s;
    private double ny;
    private double tasa_descuento_anual;
    private double CE;
    private double H;
    /*
     * variables
     */
    private double do1;
    private double L;
    private int Np;
    private double Ds;
    private double B;
    private int patron;

    private double U;
    private int Nt;
    public Intercambiadores() {
    }

    public Intercambiadores(double fluido_masico_coraza, double fluido_masico_tubo, double calor_especifico_caliente, double calor_especifico_frio, double temperatura_entrada_caliente, double temperatura_salida_caliente, double temperatura_entrada_frio, double temperatura_salda_frio, double resistencia_ensuciamiento_tubo, double resistencia_ensuciamiento_coraza, double conductividad_termica_tubo, double conductividad_termica_coraza, double densidad_tubo, double densidad_coraza, double viscosidad_dinamica_tubo, double viscosidad_dinamica_coraza, double numero_annos_vida, double tasa_descuento_anual, double costo_energia_electrica, double horas_trabajo_anno, int U) {
        this.m_s_h = fluido_masico_coraza;
        this.m_t_c = fluido_masico_tubo;
        this.cp_s_h = calor_especifico_caliente;
        this.cp_t_c = calor_especifico_frio;
        this.Ti_s_h = temperatura_entrada_caliente;
        this.To_s_h = temperatura_salida_caliente;
        this.Ti_t_c = temperatura_entrada_frio;
        this.To_t_c = temperatura_salda_frio;
        this.Rensuc_t = resistencia_ensuciamiento_tubo;
        this.Rensuc_s = resistencia_ensuciamiento_coraza;
        this.LANDA_t = conductividad_termica_tubo;
        this.LANDA_s = conductividad_termica_coraza;
        this.RO_t = densidad_tubo;
        this.RO_s = densidad_coraza;
        this.NIU_t = viscosidad_dinamica_tubo;
        this.NIU_s = viscosidad_dinamica_coraza;
        this.ny = numero_annos_vida;
        this.tasa_descuento_anual = tasa_descuento_anual;
        this.CE = costo_energia_electrica;
        this.H = horas_trabajo_anno;

        this.U = U;
        Np = 2;
        
    }

    public double getCalor_especifico_caliente() {
        return cp_s_h;
    }

    public void setCalor_especifico_caliente(double calor_especifico_caliente) {
        this.cp_s_h = calor_especifico_caliente;
    }

    public double getCalor_especifico_frio() {
        return cp_t_c;
    }

    public void setCalor_especifico_frio(double calor_especifico_frio) {
        this.cp_t_c = calor_especifico_frio;
    }

    public double getConductividad_termica_coraza() {
        return LANDA_s;
    }

    public void setConductividad_termica_coraza(double conductividad_termica_coraza) {
        this.LANDA_s = conductividad_termica_coraza;
    }

    public double getConductividad_termica_tubo() {
        return LANDA_t;
    }

    public void setConductividad_termica_tubo(double conductividad_termica_tubo) {
        this.LANDA_t = conductividad_termica_tubo;
    }

    public int getCorte_deflector() {
        return patron;
    }

    public void setCorte_deflector(int corte_deflector) {
        this.patron = corte_deflector;
    }

    public double getCosto_energia_electrica() {
        return CE;
    }

    public void setCosto_energia_electrica(double costo_energia_electrica) {
        this.CE = costo_energia_electrica;
    }

    public double getDensidad_coraza() {
        return RO_s;
    }

    public void setDensidad_coraza(double densidad_coraza) {
        this.RO_s = densidad_coraza;
    }

    public double getDensidad_tubo() {
        return RO_t;
    }

    public void setDensidad_tubo(double densidad_tubo) {
        this.RO_t = densidad_tubo;
    }

    public double getDiametro_exterior_tubo() {
        return do1;
    }

    public void setDiametro_exterior_tubo(double diametro_exterior_tubo) {
        this.do1 = diametro_exterior_tubo;
    }

    public double getDiametro_interior_coraza() {
        return Ds;
    }

    public void setDiametro_interior_coraza(double diametro_interior_coraza) {
        this.Ds = diametro_interior_coraza;
    }

    public double getEspaciamiento_baffle() {
        return B;
    }

    public void setEspaciamiento_baffle(double espaciamiento_baffle) {
        this.B = espaciamiento_baffle;
    }

    public double getFluido_masico_coraza() {
        return m_s_h;
    }

    public void setFluido_masico_coraza(double fluido_masico_coraza) {
        this.m_s_h = fluido_masico_coraza;
    }

    public double getFluido_masico_tubo() {
        return m_t_c;
    }

    public void setFluido_masico_tubo(double fluido_masico_tubo) {
        this.m_t_c = fluido_masico_tubo;
    }

    public double getHoras_trabajo_anno() {
        return H;
    }

    public void setHoras_trabajo_anno(double horas_trabajo_anno) {
        this.H = horas_trabajo_anno;
    }

    public double getLongitud_tubo() {
        return L;
    }

    public void setLongitud_tubo(double longitud_tubo) {
        this.L = longitud_tubo;
    }

    public double getNumero_annos_vida() {
        return ny;
    }

    public void setNumero_annos_vida(double numero_annos_vida) {
        this.ny = numero_annos_vida;
    }

    public double getNumero_pases_tubo() {
        return Np;
    }

    public void setNumero_pases_tubo(int numero_pases_tubo) {
        this.Np = numero_pases_tubo;
    }

    public double getResistencia_ensuciamiento_coraza() {
        return Rensuc_s;
    }

    public void setResistencia_ensuciamiento_coraza(double resistencia_ensuciamiento_coraza) {
        this.Rensuc_s = resistencia_ensuciamiento_coraza;
    }

    public double getResistencia_ensuciamiento_tubo() {
        return Rensuc_t;
    }

    public void setResistencia_ensuciamiento_tubo(double resistencia_ensuciamiento_tubo) {
        this.Rensuc_t = resistencia_ensuciamiento_tubo;
    }

    public double getTasa_descuento_anual() {
        return tasa_descuento_anual;
    }

    public void setTasa_descuento_anual(double tasa_descuento_anual) {
        this.tasa_descuento_anual = tasa_descuento_anual;
    }

    public double getTemperatura_entrada_caliente() {
        return Ti_s_h;
    }

    public void setTemperatura_entrada_caliente(double temperatura_entrada_caliente) {
        this.Ti_s_h = temperatura_entrada_caliente;
    }

    public double getTemperatura_entrada_frio() {
        return Ti_t_c;
    }

    public void setTemperatura_entrada_frio(double temperatura_entrada_frio) {
        this.Ti_t_c = temperatura_entrada_frio;
    }

    public double getTemperatura_salda_frio() {
        return To_t_c;
    }

    public void setTemperatura_salda_frio(double temperatura_salda_frio) {
        this.To_t_c = temperatura_salda_frio;
    }

    public double getTemperatura_salida_caliente() {
        return To_s_h;
    }

    public void setTemperatura_salida_caliente(double temperatura_salida_caliente) {
        this.To_s_h = temperatura_salida_caliente;
    }

    public double getViscosidad_dinamica_coraza() {
        return NIU_s;
    }

    public void setViscosidad_dinamica_coraza(double viscosidad_dinamica_coraza) {
        this.NIU_s = viscosidad_dinamica_coraza;
    }

    public double getViscosidad_dinamica_tubo() {
        return NIU_t;
    }

    public void setViscosidad_dinamica_tubo(double viscosidad_dinamica_tubo) {
        this.NIU_t = viscosidad_dinamica_tubo;
    }

    public double evaluar(double[] var) {
        do1 = var[0];
        Nt = (int) var[1];
        Ds = var[2];
        B = var[3];
        patron = (int) var[4];





        double Ctotal = 0;
        double Cinversion;
        double Coperacion;
        double Co;
        double A;
        double Q;

        double LMTD;
        double F;
        double R;
        double P;
        double di;

        double h_t;
        double Re_t;
        double Pr_t;
        double factor_friccion_darcy;
        double v_t;
        double Gs_t;
        double Astt;
        double At_t;
//        int Nt;

        double h_s;
        double De;
        double Pt;
        double Re_s;
        double Pr_s;
        double v_s;
        double Gs_s;
        double At_s;

        double petencia_bombeo;
        double deltaP_t;
        double factor_friccion_tubo;
        double deltaP_s;
        double factor_friccion_coraza;




        di = diametro_interior_tubo();

        Pt = paso_entre_tubo();

        R = coeficiente_correccion();

        P = eficiencia();

        F = factor_correccion(R, P);

        LMTD = diferencia_media_logaritmica();

        Q = carga_termica();

        A = area_superficial(Q, U, F, LMTD);

        Nt = numero_tubo();


        L = A/(Math.PI*do1*Nt);


        /*
         * coraza
         */
        At_s = area_transversal_coraza(Pt);

        Gs_s = Gs_coraza(At_s);

        v_s = velocidad_fluido_coraza(Gs_s);

        Pr_s = numero_prandtl_coraza();


        De = diametro_hidraulico_coraza(Pt);

        Re_s = numero_reynolds_coraza(v_s, De);

        h_s = coeficiente_transferencia_coraza(Re_s, Pr_s, De);

        /*
         * tubo
         */
        Pr_t = numero_prandtl_tubo();

        Astt = Astt(di);

//        Nt = numero_tubo();

        At_t = area_transversal_tubo(Astt, Nt);

        Gs_t = Gs_tubo(At_t);

        v_t = velocidad_fluido_tubo(Gs_t);

        Re_t = numero_reynolds_tubo(v_t, di);

        factor_friccion_darcy = factor_friccion_darcy(Re_t);

        h_t = coeficiente_transferencia_tubo(Re_t, di, factor_friccion_darcy, Pr_t);


        /*
         * capital inversion
         */
        U = coeficiente_global_transferencia(h_t, h_s, di);

        R = coeficiente_correccion();

        P = eficiencia();

        F = factor_correccion(R, P);

        LMTD = diferencia_media_logaritmica();

        Q = carga_termica();

        A = area_superficial(Q, U, F, LMTD);

        L = A/(Math.PI*do1*Nt);

        Cinversion = capital_inversion(A);

        /*
         * costo operacion
         */
        factor_friccion_coraza = factor_friccion_coraza(Re_s);

        deltaP_s = caida_presion_coraza(v_s, factor_friccion_coraza, De);

        factor_friccion_tubo = factor_friccion_tubo(Re_t);


//        factor_friccion_tubo = factor_friccion_darcy(Re_t);

        deltaP_t = caida_presion_tubo(v_t, di, factor_friccion_tubo);

        petencia_bombeo = petencia_bombeo(deltaP_t, deltaP_s);

        Co = Co(petencia_bombeo);

        Coperacion = costo_anual_operacion(Co);


        /*
         * general
         */
        Ctotal = Cinversion + Coperacion;
        
        Ctotal = Cinversion + Coperacion;
        if(L < 0.5){
            Ctotal += 10000;
        }

        return Ctotal;
    }

    private double diametro_interior_tubo() {
        return 0.8 * do1;
    }

    private double numero_prandtl_tubo() {
        double result = (cp_t_c * NIU_t) / LANDA_t;
        result = result * Math.pow(10, 3);
        return result;
    }

    private double Astt(double di) {
        double result = (Math.PI / 4) * Math.pow(di, 2);
        return result;
    }

    private double area_transversal_tubo(double Astt, int Nt) {
        double result = Astt * (Nt / Np);
        return result;
    }

    private double Gs_tubo(double At_t) {
        double result = m_t_c / At_t;
        return result;
    }

    private double numero_reynolds_tubo(double v_t, double di) {
        double result = (RO_t * v_t * di) / NIU_t;
        return result;
    }

    private double velocidad_fluido_tubo(double Gs_t) {
        double result = Gs_t / RO_t;
        return result;
    }

    private double factor_friccion_darcy(double Re_t) {
        double result = Math.pow((1.82 * (Math.log10(Re_t)) - 1.64), -2);

        return result;
    }

    private double coeficiente_transferencia_tubo(double Re_t, double di, double factor_friccion_darcy, double Pr_t) {
        double result;

        if (Re_t < 2300) {

            result = (LANDA_t / di) * (3.657 + (0.0677 * (Math.pow((Re_t * Pr_t * (di / L)), 1.33))) / (1 + 0.1 * Pr_t * (Math.pow((Re_t * (di / L)), 0.3))));

        } else if (Re_t >= 2300 && Re_t < 10000) {

            result = (LANDA_t / di) * ((((factor_friccion_darcy / 8) * (Re_t - 1000) * Pr_t) / (1 + 12.7 * (Math.sqrt(factor_friccion_darcy / 8)) * (Math.pow(Pr_t, 0.67) - 1))) * (1 + Math.pow((di / L), 0.67)));

        } else {

            result = (LANDA_t / do1) * 0.027 * (Math.pow(Re_t, 0.8)) * (Math.pow(Pr_t, 0.33)); // * (Math.pow((NIU_t / NIU_t), 0.14));

        }
        return result;
    }

    private double coeficiente_global_transferencia(double h_t, double h_s, double di) {
        double result = 1 / (((1 / h_s) + Rensuc_s) + ((do1 / di) * (Rensuc_t + (1 / h_t))));
        return result;
    }

    private double coeficiente_correccion() {
        double result = (Ti_s_h - To_s_h) / (To_t_c - Ti_t_c);
        return result;
    }

    private double eficiencia() {
        double result;
        result = To_t_c - Ti_t_c;
        result = result / (Ti_s_h - Ti_t_c);
        return result;
    }

    private double factor_correccion(double R, double P) {
        double result = ((Math.sqrt((Math.pow(R, 2) + 1))) / (R - 1)) * ((Math.log((1 - P) / (1 - P * R))) / (Math.log((2 - P * (R + 1 - Math.sqrt((Math.pow(R, 2) + 1)))) / (2 - P * (R + 1 + Math.sqrt((Math.pow(R, 2) + 1)))))));
        return result;
    }

    private double diferencia_media_logaritmica() {
        double result;
        double dif1 = Ti_s_h - To_t_c;
        double dif2 = To_s_h - Ti_t_c;
        result = (dif1 - dif2) / Math.log(dif1 / dif2);
        return result;
    }

    private double carga_termica() {
        double result;
        result = m_s_h * cp_s_h * (Ti_s_h - To_s_h);
        return result;
    }

    private double area_superficial(double Q, double U, double F, double LMTD) {
        double result;
        result = Q / (U * F * LMTD);
        result = result * Math.pow(10, 3);
        return result;
    }

    private double capital_inversion(double A) {
        double result;
        result = a1 + a2 * Math.pow(A, a3);
        return result;
    }

    private double paso_entre_tubo() {
        double result;
        result = 1.25 * do1;
        return result;
    }

    private double area_transversal_coraza(double Pt) {
        double result;
        result = Ds * B * (1 - (do1 / Pt));
        return result;
    }

    private double Gs_coraza(double At_s) {
        double result;
        result = m_s_h / At_s;
        return result;
    }

    private double velocidad_fluido_coraza(double Gs_s) {
        double result;
        result = Gs_s / RO_s;
        return result;
    }

    private double numero_prandtl_coraza() {
        double result;
        result = (cp_s_h * NIU_s) / LANDA_s;
        result = result * Math.pow(10, 3);
        return result;
    }

    private double diametro_hidraulico_coraza(double Np) {
        double result;
        if (patron == 45) {
            result = (1.27 / do1) * (Math.pow(Np, 2) - (0.785 * (Math.pow(do1, 2))));
        } else {
            result = (1.1 / do1) * (Math.pow(Np, 2) - (0.917 * (Math.pow(do1, 2))));
        }
        return result;
    }

    private double numero_reynolds_coraza(double v_s, double De) {
        double result;
        result = (RO_s * v_s * De) / NIU_s;
        return result;
    }

    private double coeficiente_transferencia_coraza(double Re_s, double Pr_s, double De) {
        double result;
        result = (0.36 * (Math.pow(Re_s, 0.55)) * (Math.pow(Pr_s, 0.33))) * (LANDA_s / De); //* (Math.pow((NIU_s / NIU_t), 0.14)))
        return result;
    }

    private double factor_friccion_coraza(double Re_s) {
        double result;
        result = 1.44 * (Math.pow(Re_s, -0.15));
        return result;
    }

    private double caida_presion_coraza(double v_s, double factor_friccion_coraza, double De) {
        double result;
        result = ((RO_s * Math.pow(v_s, 2)) / 2) * (factor_friccion_coraza) * (L / B) * (Ds / De);
        return result;
    }

    private double factor_friccion_tubo(double numero_reynolds_tubo) {
        double result;
        result = 0.00128 + 0.1143 * (Math.pow(numero_reynolds_tubo, -0.311));
        return result;
    }

    private double caida_presion_tubo(double v_t, double di, double factor_friccion_tubo) {
        double result;
        result = ((RO_t * Math.pow(v_t, 2)) / 2) * ((L / di) * factor_friccion_tubo + contante_p) * Np;
        return result;
    }

    private double petencia_bombeo(double delatP_t, double deltaP_s) {
        double result;
        result = (1 / eficiencia_bombeo) * ((m_t_c / RO_t) * delatP_t + (m_s_h / RO_s) * deltaP_s);
        return result;
    }

    private double Co(double petencia_bombeo) {
        double result;
        result = petencia_bombeo * CE * H;
        result = result * Math.pow(10, -3); //de llevar de w a kw la potencia de bombeo (Ce V/kw h)
        return result;
    }

    private double costo_anual_operacion(double Co) {
        double result = 0;
        for (int i = 1; i <= ny; i++) {
            result += Co / (Math.pow((1 + tasa_descuento_anual), i));
        }
        return result;
    }

    public String[] imprimir(double[] var) {
        do1 = var[0];
        Nt = (int) var[1];
        Ds = var[2];
        B = var[3];
        patron = (int) var[4];




        double Ctotal = 0;
        double Cinversion = 0;
        double Coperacion = 0;
        double Co = 0;
        double A = 0;
        double Q = 0;

        double LMTD = 0;
        double F = 0;
        double R = 0;
        double P = 0;
        double di = 0;

        double h_t = 0;
        double Re_t = 0;
        double Pr_t = 0;
        double factor_friccion_darcy = 0;
        double v_t = 0;
        double Gs_t = 0;
        double Astt = 0;
        double At_t = 0;
//        int Nt = 0;

        double h_s = 0;
        double De = 0;
        double Pt = 0;
        double Re_s = 0;
        double Pr_s = 0;
        double v_s = 0;
        double Gs_s = 0;
        double At_s = 0;

        double petencia_bombeo = 0;
        double deltaP_t = 0;
        double factor_friccion_tubo = 0;
        double deltaP_s = 0;
        double factor_friccion_coraza = 0;



        di = diametro_interior_tubo();

        Pt = paso_entre_tubo();

        R = coeficiente_correccion();

        P = eficiencia();

        F = factor_correccion(R, P);

        LMTD = diferencia_media_logaritmica();

        Q = carga_termica();

        A = area_superficial(Q, U, F, LMTD);



        Nt = numero_tubo();


        L = A/(Math.PI*do1*Nt);


        /*
         * coraza
         */
        At_s = area_transversal_coraza(Pt);

        Gs_s = Gs_coraza(At_s);

        v_s = velocidad_fluido_coraza(Gs_s);

        Pr_s = numero_prandtl_coraza();


        De = diametro_hidraulico_coraza(Pt);

        Re_s = numero_reynolds_coraza(v_s, De);

        h_s = coeficiente_transferencia_coraza(Re_s, Pr_s, De);

        /*
         * tubo
         */
        Pr_t = numero_prandtl_tubo();

        Astt = Astt(di);

//        Nt = numero_tubo();

        At_t = area_transversal_tubo(Astt, Nt);

        Gs_t = Gs_tubo(At_t);

        v_t = velocidad_fluido_tubo(Gs_t);

        Re_t = numero_reynolds_tubo(v_t, di);

        factor_friccion_darcy = factor_friccion_darcy(Re_t);



        h_t = coeficiente_transferencia_tubo(Re_t, di, factor_friccion_darcy, Pr_t);


        /*
         * capital inversion
         */
        R = coeficiente_correccion();

        P = eficiencia();

        F = factor_correccion(R, P);

        LMTD = diferencia_media_logaritmica();

        Q = carga_termica();

        U = coeficiente_global_transferencia(h_t, h_s, di);



        A = area_superficial(Q, U, F, LMTD);


        L = A/(Math.PI*do1*Nt);

        Cinversion = capital_inversion(A);

        /*
         * costo operacion
         */
        factor_friccion_coraza = factor_friccion_coraza(Re_s);

        deltaP_s = caida_presion_coraza(v_s, factor_friccion_coraza, De);

        factor_friccion_tubo = factor_friccion_tubo(Re_t);

        deltaP_t = caida_presion_tubo(v_t, di, factor_friccion_tubo);

        petencia_bombeo = petencia_bombeo(deltaP_t, deltaP_s);

        Co = Co(petencia_bombeo);

        Coperacion = costo_anual_operacion(Co);

        /*
         * general
         */
        Ctotal = Cinversion + Coperacion;
        if(L < 0.5){
            Ctotal += 10000;
        }
            
        

        String[] cadena = new String[39];
        cadena[0] = L + "";
        cadena[1] = do1 + "";
        cadena[2] = di + "";
        cadena[3] = B + "";
        cadena[4] = Ds + "";
        cadena[5] = Pt + "";
        cadena[6] = Nt + "";
        cadena[7] = Astt + "";
        cadena[8] = Gs_t + "";
        cadena[9] = At_t + "";
        cadena[10] = v_t + "";
        cadena[11] = Re_t + "";
        cadena[12] = Pr_t + "";
        cadena[13] = h_t + "";
        cadena[14] = factor_friccion_tubo + "";
        cadena[15] = deltaP_t + "";
        cadena[16] = Gs_s + "";
        cadena[17] = At_s + "";
        cadena[18] = De + "";
        cadena[19] = v_s + "";
        cadena[20] = Re_s + "";
        cadena[21] = Pr_s + "";
        cadena[22] = h_s + "";
        cadena[23] = factor_friccion_coraza + "";
        cadena[24] = deltaP_s + "";
        cadena[25] = U + "";
        cadena[26] = Q + "";
        cadena[27] = F + "";
        cadena[28] = R + "";
        cadena[29] = eficiencia_bombeo + "";
        cadena[30] = LMTD + "";
        cadena[31] = A + "";
        cadena[32] = P + "";
        cadena[33] = Cinversion + "";
        cadena[34] = Co + "";
        cadena[35] = Coperacion + "";
        cadena[36] = Ctotal + "";

        cadena[37] = patron + "";
        cadena[38] = B + "";


        


        return cadena;

    }

    public int numero_tubo() {
        int result = 0;
        double c = 0;
        double n1 = 0;

        if (patron == 30) {
            switch (Np) {
                case 1:
                    c = 0.319;
                    n1 = 2.142;
                    break;
                case 2:
                    c = 0.249;
                    n1 = 2.207;
                    break;
                case 4:
                    c = 0.175;
                    n1 = 2.285;
                    break;
                case 6:
                    c = 0.0743;
                    n1 = 2.499;
                    break;
                case 8:
                    c = 0.0365;
                    n1 = 2.675;
                    break;
            }



        } else {
            switch (Np) {
                case 1:
                    c = 0.215;
                    n1 = 2.207;
                    break;
                case 2:
                    c = 0.156;
                    n1 = 2.291;
                    break;
                case 4:
                    c = 0.158;
                    n1 = 2.263;
                    break;
                case 6:
                    c = 0.0402;
                    n1 = 2.617;
                    break;
                case 8:
                    c = 0.0331;
                    n1 = 2.643;
                    break;
            }
        }

        result = (int) (c * (Math.pow(Ds / do1, n1)));


        return result;
    }
}
