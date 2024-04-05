import csv
import pandas as pd

all_file = [["11", "12", "13", "21", "22", "23", "31", "32", "33", "41", "42", "43", "51", "52", "53", "61", "62", "63", "71", "72", "73"]]
all_pass = [["p1"]]
all_df = pd.DataFrame()

for a_pass in all_pass[0]:
    for a_file in all_file[0]:
        file_new = "csv_output/opticsDS_" + str(a_pass) + "/C12_opticsDS_" + str(a_pass)+ "_" + str(a_file) + ".csv"
        df_new=pd.read_csv(file_new)
        if not df_new.empty:
            all_df = pd.concat([all_df,df_new],axis=0,ignore_index='True') 
            #Ignore index basically makes the index numbers continuous when we are merging different files together

header=["tg_th", "tg_ph", "tg_vz", "tg_p", "gem1_r","gem1_rp","gem1_ph","gem1_php","gem1_ph_local", "sieve_r" ,"sieve_ph","rate"]
all_df.to_csv("csv_output/opticsDS_p1/C12_opticsDS_p1_merged.csv",columns=header,sep="\t")

