import math
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from polygon_selector_demo import SelectFromCollection
from scipy.stats import gaussian_kde
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model

class OPTICS:

    def __init__(self):
        self.orig = pd.DataFrame()        # data before cut
        self.sec1 = pd.DataFrame()        # sector1 data before cut
        self.sec2 = pd.DataFrame()        # sector2 data before cut
        self.sec3 = pd.DataFrame()        # sector3 data before cut
        self.sec4 = pd.DataFrame()        # sector4 data before cut
        self.sec5 = pd.DataFrame()        # sector5 data before cut
        self.sec6 = pd.DataFrame()        # sector6 data before cut
        self.sec7 = pd.DataFrame()        # sector7 data before cut

        self.selected = pd.DataFrame()    # data of selected holes

    def GenNumpyArray(self,filename):
  
        file = uproot.open(filename)
        T=file["newT"]

        geo = T.arrays(["gem1_x", "gem1_y","gem1_r","gem1_ph","gem1_px","gem1_py","gem1_pz","tg_th","tg_ph","tg_vz","tg_p","rate","sieve_r"],library="pd")  # panda dictionary
        geo = geo.loc[geo["gem1_r"]>300]

        geo["tg_ph"] = [i+2*math.pi if i<0 else i for i in geo.tg_ph]
        geo["gem1_ph"] = [i+2*math.pi if i<0 else i for i in geo.gem1_ph] 

        self.orig=geo

    def DefineSectors(self):
        rot_angle=0        # put the GEM rotation angle here

        angle_lo=[]
        angle_up=[]

        for i in range(7):
            angle_lo.append(i*2*math.pi/7 +rot_angle)
            angle_up.append((i+1)*2*math.pi/7 +rot_angle)

        print(angle_lo)
        print(angle_up)
        geo=self.orig

        self.sec1=geo.loc[(geo["gem1_ph"]<=angle_up[6]) & (geo["gem1_ph"]>=angle_lo[6])]
        self.sec2=geo.loc[(geo["gem1_ph"]<=angle_up[5]) & (geo["gem1_ph"]>=angle_lo[5])]
        self.sec3=geo.loc[(geo["gem1_ph"]<=angle_up[4]) & (geo["gem1_ph"]>=angle_lo[4])]
        self.sec4=geo.loc[(geo["gem1_ph"]<=angle_up[3]) & (geo["gem1_ph"]>=angle_lo[3])]
        self.sec5=geo.loc[(geo["gem1_ph"]<=angle_up[2]) & (geo["gem1_ph"]>=angle_lo[2])]
        self.sec6=geo.loc[(geo["gem1_ph"]<=angle_up[1]) & (geo["gem1_ph"]>=angle_lo[1])]
        self.sec7=geo.loc[(geo["gem1_ph"]<=angle_up[0]) & (geo["gem1_ph"]>=angle_lo[0])]

    def DrawHistAllSectors(self):
        fig, axs = plt.subplots(2, 4, figsize=(20,20))

        axs[0,0].hist2d(self.orig.gem1_r, self.orig.gem1_ph,(200,200),cmap=plt.cm.jet, cmin=1)
		
        axs[0,1].hist2d(self.sec1.gem1_r,self.sec1.gem1_ph,(200,200),cmap=plt.cm.jet, cmin=1)
        axs[0,1].set_title('sec1')
        axs[0,2].hist2d(self.sec2.gem1_r,self.sec2.gem1_ph,(200,200),cmap=plt.cm.jet, cmin=1)
        axs[0,2].set_title('sec2')
        axs[0,3].hist2d(self.sec3.gem1_r,self.sec3.gem1_ph,(200,200),cmap=plt.cm.jet, cmin=1)
        axs[0,3].set_title('sec3')
        axs[1,0].hist2d(self.sec4.gem1_r,self.sec4.gem1_ph,(200,200),cmap=plt.cm.jet, cmin=1)
        axs[1,0].set_title('sec4')
        axs[1,1].hist2d(self.sec5.gem1_r,self.sec5.gem1_ph,(200,200),cmap=plt.cm.jet, cmin=1)
        axs[1,1].set_title('sec5')
        axs[1,2].hist2d(self.sec6.gem1_r,self.sec6.gem1_ph,(200,200),cmap=plt.cm.jet, cmin=1)
        axs[1,2].set_title('sec6')
        axs[1,3].hist2d(self.sec7.gem1_r,self.sec7.gem1_ph,(200,200),cmap=plt.cm.jet, cmin=1)
        axs[1,3].set_title('sec7')

        plt.show()



    def SelectOneHole(self, df):

        fig, ax = plt.subplots(figsize=(10,7))

        pts=ax.scatter(df.gem1_r,df.gem1_ph)

        y_max=df.gem1_ph.max()
        y_min=df.gem1_ph.min()
        dy = (y_max-y_min)*0.1

        x_max=df.gem1_r.max()
        x_min=df.gem1_r.min()
        dx = (x_max-x_min)*0.1

        ax.set_ylim(y_min-dy, y_max+dy)
        ax.set_xlim(x_min-dx, x_max+dx)


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
        df["gem1_rp"] = (df.gem1_x*df.gem1_px+df.gem1_y*df.gem1_py)/(df.gem1_r*df.gem1_pz)     
        df["gem1_php"] = (-df.gem1_y*df.gem1_px+df.gem1_x*df.gem1_py)/(df.gem1_r*df.gem1_pz)    

        header=["tg_th","tg_ph","tg_vz","tg_p","gem1_r","gem1_rp","gem1_ph","gem1_php","sieve_r"]
        df.to_csv(filename,columns=header)

    def PolynomialRegression(self, X, y, degree, variable):  

        X_train,X_test,y_train,y_test=train_test_split(X, y, test_size=0.33, random_state=42)


        poly = PolynomialFeatures(degree)
        X_train_new=poly.fit_transform(X_train)
        X_test_new=poly.fit_transform(X_test)
 

        regression = linear_model.LinearRegression()
        model = regression.fit(X_train_new, y_train)
        y_pred = regression.predict(X_test_new)
 
        params = model.coef_
        intercept = model.intercept_

        print(intercept, params)
        score = model.score(X_test_new, y_test)
        print("score:  ", score)

        if variable == "theta":
         fig, ax = plt.subplots(1,2)
         ax[0].scatter(X_test[:,[0]],y_test[:,[0]], cmap='Greens')
         ax[0].scatter(X_test[:,[0]],y_pred[:,[0]], cmap='Reds')
         ax[0].set_ylabel('target_theta (rad)')
         ax[0].set_xlabel('gem_r(mm)')
         ax[1].scatter( X_test[:,[1]], y_test[:,[0]], cmap='Greens') 
         ax[1].scatter( X_test[:,[1]], y_pred[:,[0]], cmap='Reds')
         ax[1].set_ylabel('target_theta(rad)')
         ax[1].set_xlabel('gem_rp')
         plt.show()

        if variable == "sieve_r":
         fig, ax = plt.subplots(1,2)
         ax[0].scatter(X_test[:,[0]],y_test[:,[0]], cmap='Greens')
         ax[0].scatter(X_test[:,[0]],y_pred[:,[0]], cmap='Reds')
         ax[0].set_ylabel('sieve_r (mm)')
         ax[0].set_xlabel('gem_r(mm)')
         ax[1].scatter( X_test[:,[1]], y_test[:,[0]], cmap='Greens') 
         ax[1].scatter( X_test[:,[1]], y_pred[:,[0]], cmap='Reds')
         ax[1].set_ylabel('sieve_r(mm)')
         ax[1].set_xlabel('gem_rp')
         plt.show()

        if variable == "phi":
         fig, ax = plt.subplots(2,2)
         ax[0,0].scatter(X_test[:,[0]],y_test[:,[0]], cmap='Greens')
         ax[0,0].scatter(X_test[:,[0]],y_pred[:,[0]], cmap='Reds')
         ax[0,0].set_ylabel('target_phi (rad)')
         ax[0,0].set_xlabel('gem_r(mm)')
         ax[0,1].scatter( X_test[:,[1]], y_test[:,[0]], cmap='Greens') 
         ax[0,1].scatter( X_test[:,[1]], y_pred[:,[0]], cmap='Reds')
         ax[0,1].set_ylabel('target_phi(rad)')
         ax[0,1].set_xlabel('gem_rp')
         ax[1,0].scatter(X_test[:,[2]],y_test[:,[0]], cmap='Greens')
         ax[1,0].scatter(X_test[:,[2]],y_pred[:,[0]], cmap='Reds')
         ax[1,0].set_ylabel('target_phi (rad)')
         ax[1,0].set_xlabel('gem_phi')
         ax[1,1].scatter( X_test[:,[3]], y_test[:,[0]], cmap='Greens') 
         ax[1,1].scatter( X_test[:,[3]], y_pred[:,[0]], cmap='Reds')
         ax[1,1].set_ylabel('target_phi(rad)')
         ax[1,1].set_xlabel('gem_phip')
         plt.show()

        if variable == "momentum":
         fig, ax = plt.subplots(1,2)
         ax[0].scatter(X_test[:,[0]],y_test[:,[0]], cmap='Greens')
         ax[0].scatter(X_test[:,[0]],y_pred[:,[0]], cmap='Reds')
         ax[0].set_ylabel('target_z (mm)')
         ax[0].set_xlabel('gem_r(mm)')
         ax[1].scatter( X_test[:,[1]], y_test[:,[0]], cmap='Greens') 
         ax[1].scatter( X_test[:,[1]], y_pred[:,[0]], cmap='Reds')
         ax[1].set_ylabel('target_z(mm)')
         ax[1].set_xlabel('gem_rp')
         plt.show()
 
if __name__=='__main__':

            
     optics=OPTICS()

     flag_csv = 1
     flag_fit = 0

     if flag_csv==1:
      config = "/work/halla/moller12gev/vdoomra/optics/remoll/optics_output/slim_output/Cfoil_optics1DS_p3.root" 
      optics.GenNumpyArray(str(config))
      optics.DefineSectors()
      optics.DrawHistAllSectors()    
      optics.SelectOneHole(optics.sec7)

      hole_id="73"
      filename="output_sieve_r/optics1DS_p3_"+hole_id + ".csv"
      print(filename)
      optics.GenCSV(hole_id, filename)
 
     if flag_fit==1:    

      fit_theta = 0
      fit_phi = 0
      fit_momentum = 0
      fit_sieve_r = 1 

      if fit_theta==1 or fit_phi==1 or fit_sieve_r==1:
       all_file=["11"]
       all_target=["p3"]
       all_df = pd.DataFrame()

       for a_pass in all_target:
        for a_file in all_file:
         file_new = "output_sieve_r/optics1DS_" + str(a_pass) + "_" + str(a_file)+ ".csv"
         print(file_new)
         df_new=pd.read_csv(file_new)
         all_df = pd.concat([all_df,df_new],axis=0)

      if fit_theta==1:
       variable = "theta"
       X=all_df.iloc[:,[5,6]].values
       y=all_df.iloc[:,[1]].values
       optics.PolynomialRegression(X, y, 2, variable)

      if fit_sieve_r==1:
       variable = "sieve_r"
       X=all_df.iloc[:,[5,6]].values
       y=all_df.iloc[:,[9]].values
       optics.PolynomialRegression(X, y, 2, variable)
 
    
      if fit_phi==1:
       variable = "phi"
       X=all_df.iloc[:,[5,6,7,8]].values
       y=all_df.iloc[:,[2]].values
       optics.PolynomialRegression(X, y, 2, variable)  

      if fit_momentum==1:
        all_file=["71"]
        all_target=["p1","p2","p3"]
        all_df = pd.DataFrame()

        for a_pass in all_target:
         for a_file in all_file:
          file_new = "output/optics1DS/optics1DS_" + str(a_pass) + "_" + str(a_file) +  ".csv"
          print(file_new)
          df_new=pd.read_csv(file_new)
          all_df = pd.concat([all_df,df_new],axis=0)

        variable = "momentum"
        X=all_df.iloc[:,[5,6]].values
        y=all_df.iloc[:,[4]].values
        optics.PolynomialRegression(X, y, 3, variable)



