#Author: Vassu Doomra, Stony Brook University

import math
import ROOT
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from colorama import init, Fore, Style
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
from polygon_selector_demo import SelectFromCollection
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model

class OPTICS:

    def __init__(self):
 
        self.secNms = ['sec1', 'sec2', 'sec3', 'sec4', 'sec5', 'sec6', 'sec7']
        self.orig = pd.DataFrame()        
        self.d = {}
        self.d = {name: pd.DataFrame for name in self.secNms}
        self.selected = pd.DataFrame()
        
        self.angle_lo=[]
        self.angle_up=[]

        init()

        for i in range(7):
            self.angle_lo.append(i*2*math.pi/7)              # The entire phi region
            self.angle_up.append((i+1)*2*math.pi/7)

    def gaussian(self, x, mu, std):
       return (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-((x - mu) ** 2) / (2 * std** 2))

    def GenNumpyArray(self,filename):
  
        file = uproot.open(filename)
        T=file["newT"]

        geo = T.arrays([ "behind_sieve_k", "gem1_x", "gem1_y","gem1_r","gem1_ph","gem1_px","gem1_py","gem1_pz", "tg_th","tg_ph","tg_vz","tg_p","rate","sieve_r","sieve_ph"],library="pd")
        geo = geo.loc[geo["gem1_r"]>500]

        geo["tg_ph"] = [i+2*math.pi if i<0 else i for i in geo.tg_ph]
        geo["gem1_ph"] = [i+2*math.pi if i<0 else i for i in geo.gem1_ph]
        geo["sieve_ph"] = [i+2*math.pi if i<0 else i for i in geo.sieve_ph]

        self.orig=geo

    def DefineSectors(self):
    
        geo=self.orig

        for i in range(7):
         self.d[self.secNms[i]] = geo.loc[(geo["gem1_ph"]<=self.angle_up[i]) & (geo["gem1_ph"]>=self.angle_lo[i])]

    def DrawHistAllSectors(self):
       
        fig, bxs = plt.subplots(2, 4, figsize=(20,20))
        fig, axs = plt.subplots(2, 4, figsize=(20,20))
        fig, cxs = plt.subplots(2, 4, figsize=(20,20))
       
        for i in range(2):
         for j in range(4):

          if i*4+j == 7:
           continue

          bxs[i,j].hist(self.d[self.secNms[i*4+j]].tg_th, 200)
          bxs[i,j].set_title(self.secNms[i*4+j])
          axs[i,j].hist2d(self.d[self.secNms[i*4+j]].gem1_r,self.d[self.secNms[i*4+j]].gem1_ph,(200,200),cmap=plt.cm.jet, cmin=1)
          axs[i,j].set_title(self.secNms[i*4+j])
          cxs[i,j].hist2d(self.d[self.secNms[i*4+j]].gem1_r,self.d[self.secNms[i*4+j]].tg_th,(200,200),cmap=plt.cm.jet, cmin=1)
          cxs[i,j].set_title(self.secNms[i*4+j])

        plt.show()

    def SelectOneHole(self, df):

        fig, ax = plt.subplots(figsize=(10,7))

        pts=ax.scatter(df.gem1_r,df.gem1_ph, s=5)

        y_max=df.gem1_ph.max()
        y_min=df.gem1_ph.min()
        dy = (y_max-y_min)*0.1

        x_max=df.gem1_r.max()
        x_min=df.gem1_r.min()
        dx = (x_max-x_min)*0.1

        ax.set_ylim(y_min-dy, y_max+dy)
        ax.set_xlim(600, x_max+dx)


        selector = SelectFromCollection(ax, pts)

        print("Select points in the figure by enclosing them within a polygon.")
        print("Press the 'esc' key to start a new polygon.")
        print("Try holding the 'shift' key to move all of the vertices.")
        print("Try holding the 'ctrl' key to move a single vertex.")

        plt.show()

        selector.disconnect()
        self.selected=df.loc[df.index[selector.ind]]

    def GenCSV(self, hole_id, filename):

        df=self.selected
        
        def local_phi_transformation(x):
          for i in range(7):
            if x > self.angle_lo[i] and x < self.angle_up[i]:
             x = x - (self.angle_lo[i] + self.angle_up[i])/2
          return x
            
        df["gem1_rp"] = (df.gem1_x*df.gem1_px+df.gem1_y*df.gem1_py)/(df.gem1_r*df.gem1_pz)
        df["gem1_php"] = (-df.gem1_y*df.gem1_px+df.gem1_x*df.gem1_py)/(df.gem1_r*df.gem1_pz)
        df["gem1_ph_local"] = df.gem1_ph.apply(local_phi_transformation)

        header=["tg_th", "tg_ph", "tg_vz", "tg_p", "behind_sieve_k", "sieve_r", "sieve_ph", "gem1_r","gem1_rp","gem1_ph","gem1_php","gem1_ph_local", "gem1_px", "gem1_py", "gem1_pz","rate"]
        df.to_csv(filename,columns=header)

    def PolynomialRegression(self, X, y, degree, variable):  

        X_train,X_test,y_train,y_test=train_test_split(X, y, test_size=0.33, random_state=42)

        poly = PolynomialFeatures(degree)
        X_train_new=poly.fit_transform(X_train)
        X_test_new=poly.fit_transform(X_test)

        regression = linear_model.LinearRegression()
        model = regression.fit(X_train_new, y_train)
        y_pred_test = regression.predict(X_test_new)

        y_res_test = y_test - y_pred_test

        par = model.coef_
        intercept = model.intercept_

        print(Fore.RED + "The Fit Variable: " + Style.RESET_ALL)
        print(Fore.GREEN + variable + Style.RESET_ALL)

        print("                            ")
        score = model.score(X_test_new, y_test)
        print("score:  ", score)
        print("                            ")
        
        total_par=len(par)+1
        parameters = np.empty(total_par)
        for i in range(len(parameters)):
           if i==0:
            parameters[i] = intercept[i]
           else:
            parameters[i] = par[0][i]

        print(Fore.RED + "The Fit Parameters are: " + Style.RESET_ALL)
        print(parameters)
        print("                            ")

        n_bootstrap = 1000  # Choose an appropriate number of bootstrap samples
        intercept_samples = []
        coeff_samples = []

        for _ in range(n_bootstrap):

            indices = np.random.choice(len(X_train_new), len(X_train_new), replace=True)
            X_resampled = X_train_new[indices]
            y_resampled = y_train[indices]

            model = linear_model.LinearRegression()
            model.fit(X_resampled, y_resampled)

            # Store the parameters
            intercept_samples.append( model.intercept_)
            coeff_samples.append( model.coef_)

        intercept_samples = np.array(intercept_samples)
        coeff_samples = np.array(coeff_samples)
        intercept_uncertainties = np.std(intercept_samples, ddof=1)
        coeff_uncertainties = np.std(coeff_samples, axis=0, ddof=1)

        print(Fore.RED + "The uncertainties in the fit parameters are: " + Style.RESET_ALL)
        print(intercept_uncertainties, coeff_uncertainties[0][1:])
        print("                            ")

        varNms = ["GEM r [mm]", "GEM rp", "GEM phi [rad]", "GEM phip"] 

        fig1, bx = plt.subplots(1,4)
        fig1.canvas.manager.set_window_title(variable)
 
        for i in range(4):

          bx[i].scatter(X_test[:,[i]],y_test, s=3)
          bx[i].scatter(X_test[:,[i]],y_pred_test, s=3)
          bx[i].set_ylabel(variable)
          bx[i].set_xlabel(varNms[i])
       
        plt.show()

        hist, bin_edges, _ = plt.hist(y_res_test, bins=500, density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        initial_guess = [np.mean(y_res_test), np.std(y_res_test)]
        params, covariance = curve_fit(self.gaussian, bin_centers, hist, p0=initial_guess)
        mu, sigma = params

        x = np.linspace(min(bin_centers), max(bin_centers), 500)
        fitted_curve = self.gaussian(x, mu, sigma)
        plt.hist(y_res_test, bins=500, density=True, alpha=0.7, label='Histogram')
        plt.plot(x, fitted_curve, 'r-', label='Fitted Gaussian')
        plt.xlabel('Residuals')
        plt.legend()
        #plt.text(0.1, 0.9, f'Mean [rad/mm] = {mu:.6f}', transform=plt.gca().transAxes, fontsize=12, color='blue')
        #plt.text(0.1, 0.85, f'St. Dev. [rad/mm] = {sigma:.6f}', transform=plt.gca().transAxes, fontsize=12, color='blue')
        plt.show()

        return parameters
 
if __name__=='__main__':

    optics=OPTICS()

    flag_csv = 0
    flag_fit = 1

    if flag_csv==1:
      config = "slim_output/slim_C12_optics1US_p3.root"
      optics.GenNumpyArray(str(config))
      optics.DefineSectors()
      #optics.DrawHistAllSectors()
      optics.SelectOneHole(optics.d['sec7'])

      hole_id="73"
      filename="csv_output/optics1US_p3/C12_optics1US_p3_"+hole_id + ".csv"
      optics.GenCSV(hole_id, filename)

    if flag_fit==1:

      fitFlag   = [1, 0, 0]
      varNms    = ['theta', 'sieve_r', 'phi']
      fitDeg    = [2, 2, 2]
      ylocation = [1, 10, 2]
      par_theta = []
      par_sieve_r = []

      for i in range(3):
        if not fitFlag[i]:
            continue
        
        all_file = [["11", "12", "13", "21", "22", "23", "31", "32", "33", "41", "42", "43", "51", "52", "53", "61", "62", "63", "71", "72", "73"],
                      ["11", "12", "13", "21", "22", "23", "31", "32", "33", "41", "42", "43", "51", "52", "53", "61", "62", "63", "71", "72", "73"],
                      ["11", "12", "13", "21", "22", "23", "31", "32", "33", "41", "42", "43", "51", "52", "53", "61", "62", "63", "71", "72", "73"]]
        all_pass = [["p3"], ["p1", "p2", "p3"], ["p1", "p2", "p3"]]
        all_target = [["optics1DS", "optics1US"], ["optics1DS"], ["optics1DS"]]
        all_df = pd.DataFrame()

        for a_pass in all_pass[i]:
         for a_file in all_file[i]:
          for a_target in all_target[i]:
            file_new = "csv_output/" + str(a_target) + "_" + str(a_pass) + "/C12_" + str(a_target) + "_" + str(a_pass)+ "_" + str(a_file) + ".csv"
            print(file_new)
            df_new=pd.read_csv(file_new)
            all_df = pd.concat([all_df,df_new],axis=0, ignore_index=True)

        if i<2:

            X=all_df.iloc[:,[8,9,12,11]].values
            y=all_df.iloc[:,[ylocation[i]]].values

            if varNms[i] == 'theta':
                par_theta = optics.PolynomialRegression(X, y, fitDeg[i], varNms[i])
            if varNms[i] == 'sieve_r':
                par_sieve_r = optics.PolynomialRegression(X, y, fitDeg[i], varNms[i])

            if len(par_theta) != 0 and len(par_sieve_r) != 0:

                all_df['theta_pred'] = par_theta[0] + par_theta[1]*all_df.gem1_r + par_theta[2]*all_df.gem1_rp + par_theta[3]*all_df.gem1_r*all_df.gem1_r + par_theta[4]*all_df.gem1_r*all_df.gem1_rp + par_theta[5]*all_df.gem1_rp*all_df.gem1_rp
                all_df['sieve_r_pred'] = par_sieve_r[0] + par_sieve_r[1]*all_df.gem1_r + par_sieve_r[2]*all_df.gem1_rp + par_sieve_r[3]*all_df.gem1_r*all_df.gem1_r + par_sieve_r[4]*all_df.gem1_r*all_df.gem1_rp + par_sieve_r[5]*all_df.gem1_rp*all_df.gem1_rp
                all_df['vz_pred'] = all_df.sieve_r_pred/np.tan(all_df.theta_pred) - 250.
                all_df['vz_true'] = np.absolute(all_df.tg_vz)
                all_df['vz_true_minus_pred'] = all_df.vz_true - all_df.vz_pred

                hist, bin_edges, _ = plt.hist(all_df['vz_true_minus_pred'], bins=200, density=True)
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                initial_guess = [np.mean(all_df['vz_true_minus_pred']), np.std(all_df['vz_true_minus_pred'])]
                params, covariance = curve_fit(optics.gaussian, bin_centers, hist, p0=initial_guess)
                mu, sigma = params

                x = np.linspace(min(bin_centers), max(bin_centers), 200)
                fitted_curve = optics.gaussian(x, mu, sigma)
                plt.hist(all_df['vz_true_minus_pred'], bins=200, density=True, alpha=0.7, label='Histogram')
                plt.plot(x, fitted_curve, 'r-', label='Fitted Gaussian')
                plt.xlabel('Residuals')
                plt.legend()
                plt.text(0.1, 0.9, f'Mean [mm] = {mu:.6f}', transform=plt.gca().transAxes, fontsize=12, color='blue')
                plt.text(0.1, 0.85, f'St. Dev. [mm] = {sigma:.6f}', transform=plt.gca().transAxes, fontsize=12, color='blue')
                plt.show()

        else:

           X=all_df.iloc[:,[7,8]].values
           y=all_df.iloc[:,[ylocation[i]]].values
           optics.PolynomialRegression(X, y, fitDeg[i], varNms[i])