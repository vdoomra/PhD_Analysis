
const int nfoils = 5;
const int npar = 6;

/*

//=======The following parameters correspond to 2 degree polynomial with r, r' and phi using non-radiative events only================//

double parameter_b1[nfoils] = { -9.23952461e-03, -1.17563683e-02, -1.25830085e-02, -1.47426467e-02, -1.87400819e-02};
double parameter_b1_err[nfoils] = { 0.00017890032212018264, 0.00011426936967113017, 0.00020745396996255122, 0.00010766486590622729, 0.00010867652481347973 };

double parameter_b2[nfoils] = { 1.53183968e-04, 1.64931123e-04, 1.70337147e-04, 1.80121396e-04, 1.94781459e-04};
double parameter_b2_err[nfoils] = { 5.80357488e-07, 4.16397706e-07, 9.23236358e-07, 3.85544779e-07, 3.34477425e-07};

double parameter_b3[nfoils] = { -1.51285504, -1.58576673, -1.62709618, -1.68313969, -1.74168055};
double parameter_b3_err[nfoils] = { 4.79333108e-03, 3.54461230e-03, 8.24663706e-03, 3.42924883e-03, 2.79630902e-03 };

double parameter_b4[nfoils] = { -1.78919818e-05, 5.37363813e-04, 6.37403600e-04, 1.18675337e-03, 1.46852618e-03 };
double parameter_b4_err[nfoils] = { 1.09939361e-04, 7.25129800e-05, 1.27398618e-04, 6.86694219e-05, 5.38802023e-05 };

double parameter_b5[nfoils] = { -5.64594849e-08, -6.91163172e-08, -7.55680940e-08, -8.44918112e-08, -8.93127731e-08 };
double parameter_b5_err[nfoils] = { 1.03306881e-09, 7.73552073e-10, 2.20308016e-09, 7.23491271e-10, 5.44631827e-10 };

double parameter_b6[nfoils] = { -3.49906098e-04, -1.61688124e-04, -4.56814091e-05, 7.07816909e-05, -3.23370229e-05 };
double parameter_b6_err[nfoils] = { 2.77517410e-05, 1.90391696e-05, 5.29788372e-05, 1.79102237e-05, 1.35378306e-05 };

double parameter_b7[nfoils] = { -2.46688243e-06, -2.40423370e-06, -3.09256947e-06, -3.90124142e-06, -3.42361406e-06 };
double parameter_b7_err[nfoils] = { 4.72279195e-07, 2.80674141e-07, 5.73258203e-07, 2.54058089e-07, 1.88984817e-07 };

double parameter_b8[nfoils] = { 1.04290472e+01, 9.52891687, 8.90864565, 8.35411951, 9.50855929 };
double parameter_b8_err[nfoils] = { 1.92141243e-01, 1.21405748e-01, 3.31480746e-01, 1.18824961e-01, 9.81827591e-02 };

double parameter_b9[nfoils] = { 3.83929614e-02, 2.71101291e-02, 3.58377267e-02, 3.94271312e-02, 2.69148940e-02 };
double parameter_b9_err[nfoils] = { 7.07567296e-03, 3.79634778e-03, 7.61967688e-03, 3.40343646e-03, 2.59336463e-03 };

double parameter_b10[nfoils] = { 3.34388647e-04, -2.82728455e-04, -8.20986648e-04, -1.35175666e-03, -2.09736786e-03  };
double parameter_b10_err[nfoils] = { 1.16680789e-04, 6.58057080e-05, 1.18258651e-04, 7.07254490e-05, 5.88498009e-05  };

*/

/*

//=========The following parameters are obtained with 2 degree polynomila in r and r' from non-radiative events=========//

double parameter_b1[nfoils] = {-9.40688243e-03, -1.16444866e-02, -1.23064372e-02, -1.40510897e-02, -1.75658784e-02};
double parameter_b1_err[nfoils] = { 0.0001776722639314271, 0.00011837349777648789, 0.0001999903593484326, 0.00010182947949487774, 0.00010624208178187323 };

double parameter_b2[nfoils] = { 1.53866647e-04, 1.64716291e-04, 1.69627051e-04, 1.78109505e-04, 1.91161982e-04 };
double parameter_b2_err[nfoils] = { 5.79756299e-07, 4.15121978e-07,8.97700461e-07, 3.53771018e-07, 3.34069419e-07 };

double parameter_b3[nfoils] = { -1.51666193e+00, -1.58678410e+00, -1.62691640e+00, -1.67948631e+00, -1.73237150e+00 };
double parameter_b3_err[nfoils] = { 5.08313725e-03, 3.41045654e-03, 8.38533231e-03, 3.14528483e-03, 2.74260751e-03 };

double parameter_b4[nfoils] = { -5.59062441e-08, -6.85815806e-08, -7.44699272e-08,  -8.12408011e-08, -8.34296703e-08 };
double parameter_b4_err[nfoils] = { 1.08893358e-09, 7.65417015e-10, 2.27506380e-09, 6.82617428e-10, 5.01982909e-10 };

double parameter_b5[nfoils] = { -3.75198990e-04, -1.76012232e-04, -7.21345069e-05, -1.03617017e-06, -1.60908726e-04 };
double parameter_b5_err[nfoils] = { 2.88373735e-05, 1.88381657e-05, 5.55659956e-05, 1.71079991e-05, 1.19405758e-05 };

double parameter_b6[nfoils] = {1.06261920e+01, 9.65863893e+00, 9.14506735e+00, 8.93259193e+00, 1.05126211e+01  };
double parameter_b6_err[nfoils] = { 1.97188025e-01, 1.20354454e-01, 3.48847150e-01, 1.13150405e-01, 8.18402913e-02 };

*/

//=========The following parameters are obtained with 2 degree polynomila in r and r' from non-radiative + radiative events=========//

double parameter_b1[nfoils] = { -1.02695777e-02, -1.30587550e-02, -1.41665065e-02, -1.45986784e-02, -1.77599651e-02 };
double parameter_b1_err[nfoils] = { 0.00037718480715150896, 0.0005281920807650411, 0.0005164328749651133, 0.0001017384480751223, 0.00022111320905466965 };

double parameter_b2[nfoils] = { 1.67183164e-04, 1.75828820e-04, 1.84684934e-04, 1.79336364e-04, 1.90739459e-04 };
double parameter_b2_err[nfoils] = { 6.83048448e-06, 3.85892413e-06, 3.89711568e-06, 3.80865747e-07, 1.19003695e-06 };

double parameter_b3[nfoils] = { -1.68542632e+00, -1.70169084e+00, -1.78406466e+00, -1.67763096e+00, -1.71868038e+00 };
double parameter_b3_err[nfoils] = { 9.07901019e-02, 3.88726394e-02, 4.02281278e-02, 3.40042778e-03, 1.05959241e-02 };

double parameter_b4[nfoils] = { -1.05326498e-07, -1.06357096e-07, -1.45010056e-07, -8.35968476e-08, -8.19066911e-08 };
double parameter_b4_err[nfoils] = { 2.13144824e-08, 1.47379295e-08, 1.93825372e-08, 7.39172538e-10, 3.64613759e-09 };

double parameter_b5[nfoils] = { 8.51154631e-04, 7.61701508e-04, 1.79830710e-03, 4.73873317e-05, -2.03262026e-04 };
double parameter_b5_err[nfoils] = { 5.44808495e-04,3.79749391e-04, 5.24610279e-04, 1.81147557e-05, 9.38172957e-05   };

double parameter_b6[nfoils] = { 3.12119202e+00, 3.66096158e+00, -3.58896187e+00, 8.55453591e+00, 1.07375079e+01 };
double parameter_b6_err[nfoils] = { 3.56260278e+00, 2.53604459e+00, 3.64654556e+00, 1.16761781e-01, 6.45616359e-01 };



double x[nfoils] = {-62.45, -30., 0, 30, 62.45};

const char* parNames[npar] = { "b1", "b2", "b3", "b4", "b5", "b6"};

TGraphErrors* graph[npar] = {NULL};

void fit_parameters(){

    TCanvas* c1 = new TCanvas();
    c1->Divide(3,2);

    graph[0] = new TGraphErrors(nfoils, x, parameter_b1, 0, parameter_b1_err);
    graph[1] = new TGraphErrors(nfoils, x, parameter_b2, 0, parameter_b2_err);
    graph[2] = new TGraphErrors(nfoils, x, parameter_b3, 0, parameter_b3_err);
    graph[3] = new TGraphErrors(nfoils, x, parameter_b4, 0, parameter_b4_err);
    graph[4] = new TGraphErrors(nfoils, x, parameter_b5, 0, parameter_b5_err);
    graph[5] = new TGraphErrors(nfoils, x, parameter_b6, 0, parameter_b6_err);
    //graph[6] = new TGraphErrors(nfoils, x, parameter_b7, 0, parameter_b7_err);
    //graph[7] = new TGraphErrors(nfoils, x, parameter_b8, 0, parameter_b8_err);
    //graph[8] = new TGraphErrors(nfoils, x, parameter_b9, 0, parameter_b9_err);
    //graph[9] = new TGraphErrors(nfoils, x, parameter_b10, 0, parameter_b10_err);

    for(int ipar = 0; ipar<npar; ipar++){

        graph[ipar]->SetMarkerStyle(20);
        graph[ipar]->SetMarkerColor(kRed);
        graph[ipar]->GetXaxis()->SetTitle("Distance wrt target center [cm]");
        graph[ipar]->GetYaxis()->SetTitle(Form("%s",parNames[ipar]));
        c1->cd(ipar+1);
        graph[ipar]->Draw("AP");

    }

}  

