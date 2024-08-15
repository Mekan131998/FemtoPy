# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 18:16:54 2023

@author: Owner
"""
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
plt.close("all")

df=pd.read_csv('XPS spectra clean.txt', sep='\t', header=None, names=["KE", "BE", "BG"])

UPS=pd.read_excel('UPS spectra.xls', names=["EF", "Y_inten"])

KE=df["KE"]
BE=df["BE"]
BG=df["BG"]
df["value"]=abs(KE-BG)
EF=UPS["EF"]
Y_inten=UPS["Y_inten"]

# Filter the DataFrame based on the interval
interval_df = df[(df['BE'] >= 1180) & (df['BE'] <= 1189)]

# Find local maxima within the interval
local_maxima = (interval_df['value'] > interval_df['value'].shift(10)) & (interval_df['value'] > interval_df['value'].shift(-10))

# Display the local maxima within the interval
print(interval_df[local_maxima])


max_data = {'max_kin_en': [47.15, 307.121, 311.94, 496.498, 521.634, 628.374, 1187.70],
        'max_Y': [15660, 296895.42, 174198.24, 69606.098, 48776.734, 24310.174, 5843.489]}

max_frame = pd.DataFrame(max_data)


plt.figure(1)
plt.plot(BE, df["value"], "b")
plt.scatter(max_frame["max_kin_en"], max_frame["max_Y"], color="red", marker='x')
for i, row in max_frame.iterrows():
    plt.annotate(f'{row["max_kin_en"]:.2f} eV',
                 xy=(row['max_kin_en'], row['max_Y']),
                 xytext=(55, 1), textcoords='offset points', ha='right', va='bottom', color="green")

# plt.annotate('307.121 eV', xy=(307.121, 296895.42), xytext=(315.121, 296900.42), color='green')
# plt.annotate('311.941', xy=(315.121, 174198.24), xytext=(315.121, 174198.24), color='green')
plt.grid()
plt.xlabel("Binding energy, eV")
plt.ylabel("Intensity (counts)")
plt.tight_layout()
plt.savefig("images//Binding energy")
plt.show()
#5.39

plt.figure(2)
plt.plot(EF, Y_inten, "b")
plt.grid()
plt.xlabel("E-EF Energy, eV")
plt.ylabel("Intensity (a.u)")
plt.scatter(4.9, 50, color="red", marker='x')
plt.annotate('4.9 eV', xy=(5, 50), xytext=(5, 50), color="green")
plt.savefig("images//Intensity")
plt.show()