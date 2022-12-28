import numpy as np
import matplotlib.pyplot as plt

f=1.0
pi2=np.pi;
pi2=pi2*pi2;

f2=f*f

t0=1.0;
tt=np.linspace(0,10,201)
t2=tt-t0
t2=t2*t2
R=(1-2*pi2*f2*t2)*np.exp(-pi2*f2*t2)

W=np.fft.fft(R)
dt=tt[1]-tt[0]
Nt=len(tt)
df=1/dt/Nt
freq=np.arange(Nt)*df
fmax=1/dt*0.5

fig=plt.figure()
ax=fig.add_subplot(211)
ax.plot(tt,R)
ax.grid(True)
ax.set_title("Ricker wavelet")

bx=fig.add_subplot(212)
bx.plot(freq,np.abs(W))
bx.grid(True)
bx.set_xlim([0,fmax])
print(fmax)

plt.show()
