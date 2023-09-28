import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



path = "ovphy_plot_SLy_EOS_00.dat"

b = pd.read_csv(path,header=None,delim_whitespace=True)

for i in b.iterrows():
    print(i[1][0])