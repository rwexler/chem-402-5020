import numpy as np
import matplotlib.pyplot as plt

import mpltex

@mpltex.rsc_decorator
def my_plot(t):
    fig, ax = plt.subplots(1)
    ax.plot(t, t, label='$t$')
    ax.plot(t, t**2, label='$t^2$')
    ax.plot(t, t**3, label='$t^3$')
    ax.plot(t, t**4, label='$t^4$')
    ax.plot(t, np.log(1+t), label='$\ln(1+t)$')
    ax.plot(t, t**(1./2), label='$t^{1/2}$')
    ax.plot(t, t**(1./3), label='$t^{1/3}$')

    ax.set_xlabel('$t$')
    ax.set_ylabel('$f(t)$')
    ax.legend(loc='best', ncol=2)
    fig.tight_layout(pad=0.1)
    fig.savefig('matplotlib-rsc.png')

t = np.arange(0, 1.0+0.01, 0.01)
my_plot(t)
plt.close('all')