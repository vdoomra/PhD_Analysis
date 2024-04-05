import pandas as pd
import tensorflow as tf
import numpy as np
import warnings
from tensorflow import keras
import tensorflow_docs as tfdocs
import tensorflow_docs.plots
from tensorflow.keras import optimizers
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Activation, Dense, Dropout
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler, PolynomialFeatures
from sklearn.model_selection import GridSearchCV
from scikeras.wrappers import KerasRegressor

def predict_in_batches(model, x_data, y_data, batch_size=32):

    residuals = []
    
    for i in range(0, len(x_data), batch_size):
        x_batch = x_data[i:i+batch_size]
        y_batch_true = y_data[i:i+batch_size]
        y_batch_pred = model.predict(x_batch).flatten()

        batch_residuals = y_batch_true - y_batch_pred
        residuals.extend(batch_residuals)

        return np.array(residuals)
    
def build_model():

        model = Sequential()

        model.add(Dense(50, activation='relu', input_shape = (x_train_scaled.shape[1],)))
        model.add(Dropout(0.2))
        model.add(Dense(50, activation='relu'))
        model.add(Dropout(0.2))
        model.add(Dense(50, activation='relu'))
        model.add(Dropout(0.2))
        model.add(Dense(1))

        learning_rate = 0.1
        opt = optimizers.Adam(learning_rate)
    
        model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mse'])

        return model


if __name__=='__main__':

    tf.random.set_seed(20)
    warnings.filterwarnings("ignore")

    all_holes = ["11", "12", "13", "21", "22", "23", "31", "32", "33", "41", "42", "43", "51", "52", "53", "61", "62", "63", "71", "72", "73"]
    all_pass = ["p1", "p2" , "p3", "p4"]
    all_target = ["opticsDS"]
    all_df = pd.DataFrame()

    for a_pass in all_pass:
        for a_file in all_holes:
            for a_target in all_target:
                file_new = "csv_output/" + str(a_target) + "_" + str(a_pass) + "_non_radiative" + "/C12_" + str(a_target) + "_" + str(a_pass)+ "_" + str(a_file) + ".csv"
                df_new=pd.read_csv(file_new)
                if not df_new.empty:
                    all_df = pd.concat([all_df,df_new],axis=0, ignore_index=True)

    x=all_df.iloc[:,[5,6]].values  #r, r', phi_local and phi_prime
    y=all_df.iloc[:,[1]].values # target_theta values
    
    x_train,x_test,y_train,y_test=train_test_split(x, y, test_size=0.30, random_state=42)

    poly_features = PolynomialFeatures(degree=2)
    x_train_poly = poly_features.fit_transform(x_train)
    x_test_poly = poly_features.transform(x_test)

    scaler = MinMaxScaler()
    x_train_scaled = scaler.fit_transform(x_train_poly)
    x_test_scaled = scaler.fit_transform(x_test_poly)

    
    model = build_model()
    model.summary()

    model.fit(x_train_scaled, y_train, epochs=10, batch_size=24)

    loss, accuracy = model.evaluate(x_test_scaled, y_test)
    print(loss, accuracy)

    residuals_array = predict_in_batches(model, x_test_scaled, y_test,batch_size=24)

    plt.hist(residuals_array, bins=50, alpha=0.7, label='Histogram')
    plt.show()

    

