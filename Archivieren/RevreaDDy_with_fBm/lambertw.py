from scipy.special import lambertw
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt



a=np.linspace(-1,6)
lambertw=lambertw(a)
plt.plot(lambertw)
