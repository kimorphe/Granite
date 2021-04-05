# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm



def blank(ndat):
	Nt=ndat
	t0=0.0
	amp=np.zeros(ndat)
	time=np.arange(ndat)
	dt=1
def set_taxis(t1,t2,ndat=0):
	if ndat != 0:
		Nt=ndat
	time=np.linspace(t1,t2,Nt)
	dt=time[1]-time[0]
	df=1/(dt*Nt)
	freq=df*np.arange(Nt)



def load(nums,tb,w_6dB,mexp,dir_name=""):
	def load1(fname):
		nhead=0;
		fp=open(fname,"r")
		#fp.readline()
		#fp.readline()

		time=[]
		amp=[]
		for row in fp:
			dat=row.strip().split(",")
			time.append(float(dat[0]))
			amp.append(float(dat[1]))

		amp=amp-np.mean(amp)
		amp=amp*2.0
		time=np.array(time)*1.e6
		amp=np.array(amp)
		#amp-=np.mean(amp)
		dt=time[1]-time[0];
		Nt=len(time);
		df=1/(dt*Nt)
		freq=df*np.arange(Nt)
		fp.close()
		return time, amp, dt

	def fft(amp,dt):

		#fmax=1/dt
		#freq=np.linspace(0,fmax,Nt)
		Amp=np.fft.fft(amp);
		#Amp[1000:]=[0]*998
		#Amp[0]=0.0
		return Amp
	def Butterworth(amp,tb,w_6dB,mexp,time):
		cmp=amp
		tt=time-tb;
		Phi=1+(tt/w_6dB)**mexp
		Cmp=cmp/Phi;
		return Cmp


	V=np.array([])
	H=np.array([])
	Nx=len(nums)
	for k in nums:
		fname="scope_"+str(k)+".csv"
		if dir_name != "":
			fname=dir_name+"/"+fname
		time,amp,dt=load1(fname)
		cmp=Butterworth(amp,tb,w_6dB,mexp,time)
		Amp=fft(cmp,dt)
		V=np.hstack([V,amp])
		H=np.hstack([H,Amp])

	Nt=len(time)
	V=np.reshape(V,[Nx,Nt])
	H=np.reshape(H,[Nx,Nt])
	tmp=np.shape(V)
	Nx=tmp[0]; Nt=tmp[1];
	#self.V=V
	#self.H=H
	#self.Nx=Nx
	#self.Nt=Nt
	t1=time[0]
	t2=time[-1]
	dt=time[1]-time[0]
	dt=dt


	x1=0;
	x2=1;
	dx=0
	if Nx>1:
		dx=(x2-x1)/(Nx-1)



	return time, V, H



def load_refwv(fname):
	fp=open(fname,"r")
	#fp.readline()
	#fp.readline()

	time=[]
	amp=[]
	for row in fp:
		dat=row.strip().split(",")
		time.append(float(dat[0]))
		amp.append(float(dat[1]))

	time=np.array(time)*1.e6
	amp=np.array(amp)

	dt=time[1]-time[0];
	Nt=len(time);
	df=1/(dt*Nt)
	freq=df*np.arange(Nt)

	fp.close()
	Q=amp
	return Q,dt,Nt,time
def fft(Q,dt,Nt,normalize=True):
	W=np.fft.fft(Q)
	#W[1000:]=[0]*998
	Wmax=np.max(np.abs(W))
	df=1./(dt*Nt);
	freq=np.arange(Nt)*df
	#if normalize:
		#self.W/=self.Wmax
	#indx=np.argmax(np.abs(self.W),1)
	#self.fmax=self.freq[indx]
	#print(freq)
	return W,freq
def get_mean(Nt,Nx,amp,tb,w_6dB,mexp,time):
	def fft(amp):
		Amp=np.fft.fft(amp);
		return Amp

	def Butterworth(amp,tb,w_6dB,mexp,time):
		cmp=amp
		tt=time-tb;
		Phi=1+(tt/w_6dB)**mexp
		Cmp=cmp/Phi;
		return Cmp
	#blank(Nt)
	#set_taxis(t1,t2,Nt)
	#bsum=np.sum(self.V,0)/self.Nx

	bsum=np.zeros(Nt)
	tmp=np.zeros(Nt)
	for k in range(Nx):
		tmp=amp[k,:]
		#tmp/=np.max(tmp)
		bsum+=tmp
	bsum=bsum/Nx
	mean_amp=bsum
	mean_amp=Butterworth(mean_amp,tb,w_6dB,mexp,time)
	mean_Amp=fft(mean_amp)
	Amp=fft(bsum)
	return mean_amp, mean_Amp ,bsum ,Amp
def set_xaxis(x1,dx,npnt=0):
	x1=x1
	dx=dx
	if npnt != 0:
		Nx=npnt
	x2=x1+dx*(Nx-1)
	xcod=np.arange(Nx)*dx+x1

def wiener(Amp):
	Amax=np.max(np.abs(Amp))
	line=Amax/10

	Bmp=abs(Amp)-line
	for i in range(len(Bmp)):
		if Bmp[i] < 0:
			Bmp[i]=0.0001
	Bmp=np.array(Bmp)
	wiener=Bmp**2/abs(Amp)**2
	'''
	fig1=plt.figure()
	ax=fig1.add_subplot(111)
	ax.plot(freq,wiener)
	ax.set_xlim([0,2])
	ax.set_yscale('log')
	ax.set_ylim([1.e-4,1])

	plt.show()
	'''
	return wiener

def transfer(freq,cmp,Cmp,wiener):
	transfer=np.array([])
	for i in range(len(Cmp)):
		transfer1=Cmp[i,:]*wiener/cmp
		transfer=np.hstack([transfer,transfer1])

	transfer=transfer.reshape(len(Cmp),len(freq))

	return transfer

def mean_transfer(freq,cmp,Cmp,wiener):
	transfer=Cmp*wiener/cmp

	#print(transfer)
	return transfer

def original_transfer(freq,cmp,Cmp):
	transfer2=Cmp/cmp

	#print(transfer)
	return transfer2


def refButterworth(reference,tb,w_6dB,mexp,time):
	tt=time-tb;
	Phi=1+(tt/w_6dB)**mexp
	reference1=reference/Phi;
	return reference1

def ref_amax(amp,time):
	amax=[]
	max_time=[]
	bmp=amp
	imax=np.argmax(-bmp)
	time_max=time[imax]
	tmax=bmp[imax]
	max_time.append(time_max)
	amax.append(tmax)

	max_time=np.array(max_time)
	amax=np.array(amax)
	return max_time, amax

def amax(amp,time,Nx):
	amax=[]
	max_time=[]
	for k in range(Nx):
		bmp=amp[k,:]
		imax=np.argmax(-bmp)
		time_max=time[imax]
		tmax=bmp[imax]
		max_time.append(time_max)
		amax.append(tmax)

	max_time=np.array(max_time)
	amax=np.array(amax)
	return max_time, amax

def max_variance(time,amp,mean_time,mean_amp,Nx):
	times=[]
	cmax=[]
	for k in range(Nx):
		bmp=amp[k]
		ttime=time[k]
		tvariance=(ttime-mean_time)**2
		times.append(tvariance)
		avariance=(bmp-mean_amp)**2
		cmax.append(avariance)

	cmax=np.array(cmax)
	times=np.array(times)
	cmax=np.std(cmax)
	times=np.std(times)
	return cmax, times

def gdelays(transfer1,freq,Nx):
	fmin=0.5
	fmax=1.6
	phis=[]
	tgs=[]
	cgs=[]
	L=3.42
	for k in range(Nx):
		phi=np.unwrap(np.angle(transfer1[k,:]))
		#dw=2.*np.pi*self.df
		nf1=np.argmin(np.abs((fmin-freq)))
		nf2=np.argmin(np.abs((fmax-freq)))
		tg,t0=np.polyfit(2.*np.pi*freq[nf1:nf2],phi[nf1:nf2],1)
		cg=L/tg
		phis.extend(phi)
		tgs.append(tg)
		cgs.append(cg)
	phis=np.array(phis)
	phis=phis.reshape(Nx,2000)
	tgs=np.array(tgs)
	cgs=np.array(cgs)
	#phis=-phis
	tgs=-tgs
	cgs=-cgs
	return phis, tgs ,cgs

def mean_gdelays(transfer,freq):
	fmin=0.5
	fmax=1.6

	phi=np.unwrap(np.angle(transfer))
	#dw=2.*np.pi*self.df
	nf1=np.argmin(np.abs((fmin-freq)))
	nf2=np.argmin(np.abs((fmax-freq)))
	tg,t0=np.polyfit(2.*np.pi*freq[nf1:nf2],phi[nf1:nf2],1)
	phi=np.array(phi)
	tg=-tg
	return phi, tg

def refft(Amp, Nx):
	reamps=np.array([])
	for i in range(Nx):
		re_amp = np.fft.ifft(Amp[i,:])
		re_amp = re_amp.real
		reamps = np.hstack([reamps,re_amp])

	reamps=reamps.reshape(Nx,len(Amp[0,:]))
	return reamps

def mean_refft(mean_Amp):
	re_amp = np.fft.ifft(mean_Amp)
	mean_f = re_amp.real
	return mean_f

def fnc_std(func,Nx):
	mayu=np.std(func,0)
	return mayu

def func_phase(Amp,Nx):
	phases = []
	for i in range(Nx):
		phase1 = [np.arctan2(float(c.imag), float(c.real)) for c in Amp[i,:]]
		phases.extend(phase1)
	phases = np.array(phases)
	#print(phases)
	#print(len(phases))
	phases = phases.reshape(Nx,len(Amp[0,:]))

	vsum=np.zeros(len(Amp[0,:]))
	for k in range(Nx):
		vmp=phases[k,:]
		vsum+=vmp
	vsum=vsum/Nx
	mean_phases=vsum
	return phases, mean_phases

def mean_func_phase(Amp):
	phase = [np.arctan2(float(c.imag), float(c.real)) for c in Amp]
	return phase

def tile(transfer,freq,Nx):
	fmin=0.5
	fmax=1.5
	transtiles=[]
	for k in range(Nx):
		nf1=np.argmin(np.abs((fmin-freq)))
		nf2=np.argmin(np.abs((fmax-freq)))
		tg,t0=np.polyfit(freq[nf1:nf2],transfer[k,nf1:nf2],1)
		transtiles.append(tg)
	transtiles=np.array(transtiles)
	#transtiles=-transtiles
	return transtiles

def transfer_max(transfer, freq, Nx):
	min=0.5
	fmax=1.5
	transmaxs=[]
	for k in range(Nx):
		nf1=np.argmin(np.abs((fmin-freq)))
		nf2=np.argmin(np.abs((fmax-freq)))
		max1=np.argmax(transfer[k,:nf2])
		transmax=transfer[k,max1]
		transmaxs.append(transmax)
	transmaxs=np.array(transmaxs)
	#transtiles=-transtiles
	return transmaxs

def mean_transfer_max(transfer, freq):
	min=0.5
	fmax=1.5
	nf2=np.argmin(np.abs((fmax-freq)))
	max1=np.argmax(transfer[:nf2])
	transmax=transfer[max1]
	#transtiles=-transtiles
	return transmax, max1

def count(transfer, Nx, nf2):
	x=np.linspace(0,1.01,101)
	transfer=transfer[:,:nf2]
	bs=[]
	for i in range(200):
		cs=[]
		for k in range(100):
			c=np.count_nonzero((x[k]<transfer[:,i])&(transfer[:,i]<x[k+1]))
			cs.append(c)
		csmax=np.max(cs)
		cs=cs/csmax
		bs.extend(cs)
	bs=np.array(bs)

	bs=bs.reshape(200,100)
	#bsmax=np.max(bs[100,:])
	#bs=bs/bsmax
	print(bs.shape)
	return bs

def count1(transfer, Nx, nf2):
	x=np.linspace(2.5,-20.1,226)
	transfer=transfer[:,:nf2]
	bs=[]
	for i in range(200):
		cs=[]
		for k in range(225):
			c=np.count_nonzero((x[k]<transfer[:,i])&(transfer[:,i]<x[k+1]))
			cs.append(c)
		csmax=np.max(cs)
		cs=cs/csmax
		bs.extend(cs)
	bs=np.array(bs)

	bs=bs.reshape(200,225)
	#bsmax=np.max(bs[100,:])
	#bs=bs/bsmax
	print(bs.shape)
	return bs

def hosei(transfer, Nx, c):
	#transfer=transfer[:,:nf2]
	bs=[]
	for i in range(Nx):
		trans=transfer[i,:]*c
		bs.extend(trans)
	bs=np.array(bs)
	bs=bs.reshape(Nx,2000)
	return bs

def gdelay1(transfer,freq,Nx,df):
	fmin=0.5
	fmax=1.5
	PI=np.pi
	phis=[]
	tgs=[]
	cgs=[]
	L=3.42
	for k in range(Nx):
		phi=np.unwrap(np.angle(transfer[k,:]))
		#dw=2.*np.pi*self.df
		nf1=np.argmin(np.abs((fmin-freq)))
		nf2=np.argmin(np.abs((fmax-freq)))
		tg=-np.diff(phi)/(2*df*PI) #	Group Delay
		for i in range(1999):
			cg=L/tg[i]
			cgs.append(cg)
		phis.extend(phi)
		tgs.extend(tg)
		#cgs.extend(cg)
	phis=np.array(phis)
	phis=phis.reshape(Nx,2000)
	tgs=np.array(tgs)
	cgs=np.array(cgs)
	#phis=-phis
	tgs=tgs.reshape(Nx,1999)
	cgs=cgs.reshape(Nx,1999)
	return  tgs ,cgs

def mean_gdelay1(transfer,freq,df):
	fmin=0.5
	fmax=1.5
	PI=np.pi
	phis=[]
	tgs=[]
	cgs=[]
	L=3.42
	phi=np.unwrap(np.angle(transfer))
	#dw=2.*np.pi*self.df
	nf1=np.argmin(np.abs((fmin-freq)))
	nf2=np.argmin(np.abs((fmax-freq)))
	tg=-np.diff(phi)/(2*df*PI) #	Group Delay
	for i in range(1999):
		cg=L/tg[i]
		cgs.append(cg)
	phis.extend(phi)
	tgs.extend(tg)
	#cgs.extend(cg)
	phis=np.array(phis)
	tgs=np.array(tgs)
	cgs=np.array(cgs)
	return  tgs ,cgs

def esalab(transfer, guntien):
	wao=[]
	for k in range(199):
		guntien10=guntien[:,k]
		guntien10=guntien1.T
		transfer10=abs(transfer[:,k+1])
		transfer10=transfer1.T
		esalab1=np.cov(guntien10,transfer10)
		wao.append(esalab1[0,1])

	wao=np.array(wao)
	print(len(wao))
	return  wao





#Main program
#------------ B-scan Data Description -------------------------------
sekinums=range(590)  #sekiei File Numbers (scope_***.csv)
karinums=range(825)  #kariumu File Numbers (scope_***.csv)
natonums=range(549)  #natoriumu File Numbers (scope_***.csv)
Nx=len(sekinums)+len(karinums)+len(natonums)    # Number of waveforms
x1=0; x2= 30; # Observation line (start,end)
dx=(x2-x1)/(Nx-1) #pitch
L=3  # travel distance [mm]
reftb=11.4 # approximate arrival time [micro sec]
tb=12.00 # approximate arrival time [micro sec]
w_6dB=1.2; # 6dB half pulse width [micro sec]
mexp=4;     # Butterworth window parameter (exponent)

seki_c=2.083
kari_c=1.887
nato_c=2.027
fmin=0.5
fmax=1.6

#------------  Reference Signal  --------------------------
fname_ref="tansyokushi1MHz/1MHznew.csv" # file name
reference,dt,Nt,time=load_refwv(fname_ref) #load reference signal
reference1=refButterworth(reference,reftb,w_6dB,mexp,time)
spectrum_reference,freq=fft(reference1,dt,Nt) # Butterworth filtering

df=1/dt/Nt
#num_dir=np.arange(1,11,1);
transfers=[]; freqs=[]
mean_transfers=[]

amp=[]
Amp=[]

#------------  Load Files  ----------------------------------
dir_name='seki'  # Data Directory
time,amp1,Amp1=load(sekinums,tb,w_6dB,mexp,dir_name=dir_name) # load waveforms and time axis
amp.extend(amp1)
Amp.extend(Amp1)

dir_name='kari'  # Data Directory
time,amp2,Amp2=load(karinums,tb,w_6dB,mexp,dir_name=dir_name) # load waveforms and time axis
amp.extend(amp2)
Amp.extend(Amp2)

dir_name='nato'  # Data Directory
time,amp3,Amp3=load(natonums,tb,w_6dB,mexp,dir_name=dir_name) # load waveforms and time axis
amp.extend(amp3)
Amp.extend(Amp3)

amp=np.array(amp)
Amp=np.array(Amp)

#print(Amp[0])
#print(len(Amp[0]))
#print(time.shape)

#------------  Draw B-scan  -----------------------------
mean_amp,mean_Amp,mean_amp_original,mean_Amp_original=get_mean(Nt,Nx,amp,tb,w_6dB,mexp,time)  # get mean waveform

mean_ampx,mean_Ampx,mean_ampx_original,mean_Ampx_original=get_mean(Nt,len(sekinums),amp1,tb,w_6dB,mexp,time)
mean_ampy,mean_Ampy,mean_ampy_original,mean_Ampy_original=get_mean(Nt,len(karinums),amp2,tb,w_6dB,mexp,time)
mean_ampz,mean_Ampz,mean_ampz_original,mean_Ampz_original=get_mean(Nt,len(natonums),amp3,tb,w_6dB,mexp,time)

#------------  Wiener Filter  --------------------------------
wiener=wiener(mean_Amp)
#mean_wiener=wiener(mean_Amp)

#------------  transfer function  ------------------------
#transfer1=transfer(freq,spectrum_reference,Amp,wiener)
transfer1=transfer(freq,spectrum_reference,Amp,wiener)
transfers.extend(transfer1)
transfers=np.array(transfers)

seki_transfer=transfers[0:len(sekinums),:]
kari_transfer=transfers[len(sekinums):len(sekinums)+len(karinums),:]
nato_transfer=transfers[len(sekinums)+len(karinums):len(sekinums)+len(karinums)+len(natonums),:]

seki_transfer1=hosei(seki_transfer, len(sekinums),seki_c)
kari_transfer1=hosei(kari_transfer, len(karinums),kari_c)
nato_transfer1=hosei(nato_transfer, len(natonums),nato_c)

#------------  transfer function std  ----------------------
#seki_std=fnc_std(seki_transfer,len(sekinums))
#kari_std=fnc_std(kari_transfer,len(karinums))
#nato_std=fnc_std(nato_transfer,len(natonums))

#------------  mean transfer function  -----------
mean_transfer1=mean_transfer(freq,spectrum_reference,mean_Amp,wiener)
mean_transfers.extend(mean_transfer1)

mean_transferx=mean_transfer(freq,spectrum_reference,mean_Ampx,wiener)
mean_transferx=mean_transferx*seki_c
mean_transfers.extend(mean_transferx)
mean_transfery=mean_transfer(freq,spectrum_reference,mean_Ampy,wiener)
mean_transfery=mean_transfery*kari_c
mean_transfers.extend(mean_transfery)
mean_transferz=mean_transfer(freq,spectrum_reference,mean_Ampz,wiener)
mean_transferz=mean_transferz*nato_c
mean_transfers.extend(mean_transferz)

#Amp[:,Nt//2:]=0*(Nt//2)
#freq[Nt//2:]=0*(Nt//2)

mean_transfers=np.array(mean_transfers)
mean_transfers=mean_transfers.reshape(4,Nt)
#print(len(transfers))

#original_seki_transfer=original_transfer(freq,spectrum_reference,mean_Ampx_original)
#original_kari_transfer=original_transfer(freq,spectrum_reference,mean_Ampy_original)
#original_nato_transfer=original_transfer(freq,spectrum_reference,mean_Ampz_original)

#------------  高速逆フーリエ変換  実部の値のみ取り出し  -----------
reamp = refft(transfers,Nx)
mean_reamp_total = mean_refft(mean_transfers[0])
mean_reamp_seki = mean_refft(mean_transfers[1])
mean_reamp_kari = mean_refft(mean_transfers[2])
mean_reamp_nato = mean_refft(mean_transfers[3])

#-------------  Peak max variance  -------------------
max_time,max_amp=amax(amp,time,Nx)
mean_max_time,mean_max_amp=ref_amax(mean_amp,time)
time_variance,max_variance=max_variance(max_time,max_amp,mean_max_time,mean_max_amp,Nx)
#------------  transfer function Declination  --------------------------
declination,guntien,gunsokudo=gdelays(transfers,freq,Nx)
seki_mean_declination,seki_mean_guntien=mean_gdelays(mean_transferx,freq)
kari_mean_declination,kari_mean_guntien=mean_gdelays(mean_transfery,freq)
nato_mean_declination,nato_mean_guntien=mean_gdelays(mean_transferz,freq)

seki_guntien=guntien[0:len(sekinums)]
kari_guntien=guntien[len(sekinums):len(sekinums)+len(karinums)]
nato_guntien=guntien[len(sekinums)+len(karinums):len(sekinums)+len(karinums)+len(natonums)]
'''
print(np.median(seki_guntien))
print(np.median(kari_guntien))
print(np.median(nato_guntien))
print(np.mean(seki_guntien))
print(np.mean(kari_guntien))
print(np.mean(nato_guntien))
print(np.std(seki_guntien))
print(np.std(kari_guntien))
print(np.std(nato_guntien))
'''
seki_declination=declination[0:len(sekinums),:]
kari_declination=declination[len(sekinums):len(sekinums)+len(karinums),:]
nato_declination=declination[len(sekinums)+len(karinums):len(sekinums)+len(karinums)+len(natonums),:]

seki_gunsokudo=gunsokudo[0:len(sekinums)]
kari_gunsokudo=gunsokudo[len(sekinums):len(sekinums)+len(karinums)]
nato_gunsokudo=gunsokudo[len(sekinums)+len(karinums):len(sekinums)+len(karinums)+len(natonums)]
'''
print(np.median(seki_gunsokudo))
print(np.median(kari_gunsokudo))
print(np.median(nato_gunsokudo))
print(np.mean(seki_gunsokudo))
print(np.mean(kari_gunsokudo))
print(np.mean(nato_gunsokudo))
print(np.std(seki_gunsokudo))
print(np.std(kari_gunsokudo))
print(np.std(nato_gunsokudo))
'''
#------------  phase function  -----------------------
#seki_phase, mean_sekis_phase = func_phase(seki_transfer,len(sekinums))
#kari_phase, mean_karis_phase = func_phase(kari_transfer,len(karinums))
#nato_phase, mean_natos_phase = func_phase(nato_transfer,len(natonums))

#mean_seki_phase = mean_func_phase(mean_transferx)
#mean_kari_phase = mean_func_phase(mean_transfery)
#mean_nato_phase = mean_func_phase(mean_transferz)

#------------  transfer probability  -----------------------
seki_transfer_tilt = tile(abs(seki_transfer), freq, len(sekinums))
kari_transfer_tilt = tile(abs(kari_transfer), freq, len(karinums))
nato_transfer_tilt = tile(abs(nato_transfer), freq, len(natonums))

#------------  transfer max  -----------------------
seki_transfermax = transfer_max(abs(seki_transfer), freq, len(sekinums))
kari_transfermax = transfer_max(abs(kari_transfer), freq, len(karinums))
nato_transfermax = transfer_max(abs(nato_transfer), freq, len(natonums))

#------------  mean transfer max  -----------------------
seki_mean_transfermax,maxx2 = mean_transfer_max(abs(mean_transferx), freq)
kari_mean_transfermax,maxy2 = mean_transfer_max(abs(mean_transfery), freq)
nato_mean_transfermax,maxz2 = mean_transfer_max(abs(mean_transferz), freq)
'''
print(seki_mean_transfermax)
print(kari_mean_transfermax)
print(nato_mean_transfermax)
print(maxx2)
print(maxy2)
print(maxz2)


print(abs(mean_transferx[150]))
print(abs(mean_transfery[150]))
print(abs(mean_transferz[150]))
'''
#------------  transfer probability + transfer max  -----------------------
z_seki = []
z_seki.extend(seki_transfer)
z_seki = np.array(z_seki)
#z_seki = z_seki.reshape(len(sekinums)+1,len(freq))
z_kari = []
z_kari.extend(kari_transfer)
z_kari = np.array(z_kari)
#z_kari = z_kari.reshape(2,len(karinums))
z_nato = []
z_nato.extend(nato_transfer)
z_nato = np.array(z_nato)
#z_nato = z_nato.reshape(2,len(natonums))

#------------  transfer probability + transfer max  -----------------------
fmin=0.5
fmax=2.0
nf1=np.argmin(np.abs((fmin-freq)))
nf2=np.argmin(np.abs((fmax-freq)))
seki_count = count(abs(seki_transfer1), len(sekinums),nf2 )
kari_count = count(abs(kari_transfer1), len(karinums),nf2 )
nato_count = count(abs(nato_transfer1), len(natonums),nf2 )
#seki_count1 = count1(seki_declination, len(sekinums),nf2 )
#kari_count1 = count1(kari_declination, len(karinums),nf2 )
#nato_count1 = count1(nato_declination, len(natonums),nf2 )
seki_count=seki_count.T
seki_count=seki_count[::-1,:]
kari_count=kari_count.T
kari_count=kari_count[::-1,:]
nato_count=nato_count.T
nato_count=nato_count[::-1,:]

#------------  transfer std  -----------------------
seki_std=np.std(abs(seki_transfer1[:,:nf2]),0)
kari_std=np.std(abs(kari_transfer1[:,:nf2]),0)
nato_std=np.std(abs(nato_transfer1[:,:nf2]),0)
a=np.argmax(seki_std)
b=np.argmax(kari_std)
c=np.argmax(nato_std)

print(a)
print(b)
print(c)
print(seki_std.max())
print(kari_std.max())
print(nato_std.max())

#------------  delay  -----------------------
freq1=freq[1:]
guntien1,gunsokudo1=gdelay1(seki_transfer1,freq,len(sekinums),df)
guntien2,gunsokudo2=gdelay1(kari_transfer1,freq,len(karinums),df)
guntien3,gunsokudo3=gdelay1(nato_transfer1,freq,len(natonums),df)
mean_guntien1,mean_gunsokudo1=mean_gdelay1(mean_transferx,freq,df)
mean_guntien2,mean_gunsokudo2=mean_gdelay1(mean_transfery,freq,df)
mean_guntien3,mean_gunsokudo3=mean_gdelay1(mean_transferz,freq,df)
guntien1=guntien1[:,50:nf2]
guntien2=guntien2[:,50:nf2]
guntien3=guntien3[:,50:nf2]
gunsokudo1=gunsokudo1[:,50:nf2]
gunsokudo2=gunsokudo2[:,50:nf2]
gunsokudo3=gunsokudo3[:,50:nf2]
flatguntien_seki=guntien1.flatten()
flatguntien_kari=guntien2.flatten()
flatguntien_nato=guntien3.flatten()
flatgunsokudo_seki=gunsokudo1.flatten()
flatgunsokudo_kari=gunsokudo2.flatten()
flatgunsokudo_nato=gunsokudo3.flatten()
print(np.median(flatgunsokudo_seki))
print(np.median(flatgunsokudo_kari))
print(np.median(flatgunsokudo_nato))
print(np.mean(flatgunsokudo_seki))
print(np.mean(flatgunsokudo_kari))
print(np.mean(flatgunsokudo_nato))
print(np.std(flatgunsokudo_seki))
print(np.std(flatgunsokudo_kari))
print(np.std(flatgunsokudo_nato))


'''
#------------  共分散  -----------------------
seki_tienbunsan=esalab(seki_transfer1, guntien1)
kari_tienbunsan=esalab(kari_transfer1, guntien2)
nato_tienbunsan=esalab(nato_transfer1, guntien3)

seki_gunbunsan=esalab(seki_transfer1, gunsokudo1)
kari_gunbunsan=esalab(kari_transfer1, gunsokudo2)
nato_gunbunsan=esalab(nato_transfer1, gunsokudo3)

#------------  show python program  -----------
freq2=freq1[0:199]
fig1=plt.figure()
ax=fig1.add_subplot(111)
ax.plot(freq2,seki_tienbunsan,'b',label="seki")
ax.plot(freq2,kari_tienbunsan,'r',label="kari")
ax.plot(freq2,nato_tienbunsan,'g',label="nato")
ax.grid(True)

fsz=10; # fontsize
ax.set_xlabel("frequency(MHz)",fontsize=fsz)
ax.set_ylabel("Distributed",fontsize=fsz)
ax.legend()
ax.tick_params(labelsize=fsz)
ax.set_xlim([20,170])
#ax.set_ylim([0,15])
fig1.savefig("bunsantien.png",bbox_inches="tight")





fig, ax = plt.subplots()
im=ax.imshow(seki_count,vmin=0,vmax=1,cmap="jet")
#x_ticks=np.arange(0,200,20)
#y_ticks=np.arange(0,150,15)

#ax.set_xticks(x_ticks)
#ax.set_yticks(y_ticks)
#heatmap = ax.pcolor(seki_count, cmap=plt.cm.Blues,vmax=1,vmin=0)
#ax.set_ylim(50,150)

fig1, ax1 = plt.subplots()
im=ax1.imshow(kari_count,vmin=0,vmax=1,cmap="jet")
#x_ticks=np.arange(0,200,20)
#y_ticks=np.arange(0,150,15)

#ax1.set_xticks(x_ticks)
#ax1.set_yticks(y_ticks)

fig2, ax2 = plt.subplots()
im=ax2.imshow(nato_count,vmin=0,vmax=1,cmap="jet")
#x_ticks=np.arange(0,200,20)
#y_ticks=np.arange(0,150,15)

#ax2.set_xticks(x_ticks)
#ax2.set_yticks(y_ticks)

plt.show()

mean_guntien4=np.mean(guntien1,0)
mean_guntien5=np.mean(guntien2,0)
mean_guntien6=np.mean(guntien3,0)
mean_gunsokudo4=np.mean(gunsokudo1,0)
mean_gunsokudo5=np.mean(gunsokudo2,0)
mean_gunsokudo6=np.mean(gunsokudo3,0)
#------------  散布図 振幅　群遅延  -----------
sekit=np.mean(abs(seki_transfer1),0)
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
for i in range(len(sekinums)-1):
	ax1.scatter(guntien1[i,49],abs(seki_transfer1[i,50]),s=10, c='b',marker='.')
ax1.scatter(guntien1[len(sekinums)-1,49],abs(seki_transfer1[len(sekinums)-1,50]),label='0.5MHz',s=10, c='b',marker='.')
for i in range(len(sekinums)-1):
	ax1.scatter(guntien1[i,99],abs(seki_transfer1[i,100]),s=10,c='red', marker='.')
ax1.scatter(guntien1[len(sekinums)-1,99],abs(seki_transfer1[len(sekinums)-1,100]),label='1.0MHz',s=10,c='red', marker='.')
for i in range(len(sekinums)-1):
	ax1.scatter(guntien1[i,149],abs(seki_transfer1[i,150]),s=10,c='g', marker='.')
ax1.scatter(guntien1[len(sekinums)-1,149],abs(seki_transfer1[len(sekinums)-1,150]),label='1.5MHz',s=10,c='g', marker='.')
#ax1.scatter(mean_guntien1[49:149],abs(mean_transferx[50:150]),c='y', label='mean transfer',marker='.')
ax1.scatter(mean_guntien4[49],sekit[50], s=200,c='y',label='mean_0.5MHz',marker='.')
ax1.scatter(mean_guntien4[99],sekit[100],s=200,c=(0.3, 0.2, 0.8) , label='mean_1.0MHz',marker='.')
ax1.scatter(mean_guntien4[149],sekit[150],s=200, c='m',label='mean_1.5MHz',marker='.')
ax1.set_title('Covariance')
ax1.set_xlabel('group_delay[μsec]')
ax1.set_ylabel('amplitude')
ax1.set_xlim([0,2])
ax1.set_ylim([0,1.0])
ax1.grid(True)
ax1.legend()
fig1.savefig("kyoubunsan_seki.png",bbox_inches="tight")

karit=np.mean(abs(kari_transfer1),0)
fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
for i in range(len(karinums)-1):
	ax2.scatter(guntien2[i,49],abs(kari_transfer1[i,50]),s=10, c='b',marker='.')
ax2.scatter(guntien2[len(karinums)-1,49],abs(kari_transfer1[len(karinums)-1,50]),label='0.5MHz',s=10, c='b',marker='.')
for i in range(len(karinums)-1):
	ax2.scatter(guntien2[i,99],abs(kari_transfer1[i,100]),s=10,c='red', marker='.')
ax2.scatter(guntien2[len(karinums)-1,99],abs(kari_transfer1[len(karinums)-1,100]),label='1.0MHz',s=10,c='red', marker='.')
for i in range(len(karinums)-1):
	ax2.scatter(guntien2[i,149],abs(kari_transfer1[i,150]),s=10,c='g', marker='.')
ax2.scatter(guntien2[len(karinums)-1,149],abs(kari_transfer1[len(karinums)-1,150]),label='1.5MHz',s=10,c='g', marker='.')
#ax2.scatter(mean_guntien2[49:149],abs(mean_transfery[50:150]),c='y',label='mean transfer', marker='.')
ax2.scatter(mean_guntien5[49],karit[50],s=200,c='y',label='mean_0.5MHz', marker='.')
ax2.scatter(mean_guntien5[99],karit[100],s=200,c=(0.3, 0.2, 0.8) ,label='mean_1.0MHz', marker='.')
ax2.scatter(mean_guntien5[149],karit[150],s=200,c='m',label='mean_1.5MHz', marker='.')
ax2.set_title('Covariance')
ax2.set_xlabel('group_delay[μsec]')
ax2.set_ylabel('amplitude')
ax2.set_xlim([0,2])
ax2.set_ylim([0,1.0])
ax2.grid(True)
ax2.legend()
fig2.savefig("kyoubunsan_kari.png",bbox_inches="tight")

natot=np.mean(abs(nato_transfer1),0)
fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
for i in range(len(natonums)-1):
	ax3.scatter(guntien3[i,49],abs(nato_transfer1[i,50]),s=10, c='b',marker='.')
ax3.scatter(guntien3[len(natonums)-1,49],abs(nato_transfer1[len(natonums)-1,50]),label='0.5MHz',s=10, c='b',marker='.')
for i in range(len(natonums)-1):
	ax3.scatter(guntien3[i,99],abs(nato_transfer1[i,100]),s=10,c='red', marker='.')
ax3.scatter(guntien3[len(natonums)-1,99],abs(nato_transfer1[len(natonums)-1,100]),label='1.0MHz',s=10,c='red', marker='.')
for i in range(len(natonums)-1):
	ax3.scatter(guntien3[i,149],abs(nato_transfer1[i,150]),s=10,c='g', marker='.')
ax3.scatter(guntien3[len(natonums)-1,149],abs(nato_transfer1[len(natonums)-1,150]),label='1.5MHz',s=10,c='g', marker='.')
#ax3.scatter(mean_guntien3[49:149],abs(mean_transferz[50:150]),c='y',label='mean transfer', marker='.')
ax3.scatter(mean_guntien6[49],natot[50],s=200,c='y',label='mean_0.5MHz', marker='.')
ax3.scatter(mean_guntien6[99],natot[100],s=200,c=(0.3, 0.2, 0.8) ,label='mean_1.0MHz', marker='.')
ax3.scatter(mean_guntien6[149],natot[150],s=200,c='m',label='mean_1.5MHz', marker='.')
ax3.set_title('Covariance')
ax3.set_xlabel('group_delay[μsec]')
ax3.set_ylabel('amplitude')
ax3.set_xlim([0,2])
ax3.set_ylim([0,1.0])
ax3.grid(True)
ax3.legend()
fig3.savefig("kyoubunsan_nato.png",bbox_inches="tight")

#郡速度

fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
for i in range(len(sekinums)-1):
	ax1.scatter(gunsokudo1[i,49],abs(seki_transfer1[i,50]),s=10, c='b',marker='.')
ax1.scatter(gunsokudo1[len(sekinums)-1,49],abs(seki_transfer1[len(sekinums)-1,50]),label='0.5MHz', s=10,c='b',marker='.')
for i in range(len(sekinums)-1):
	ax1.scatter(gunsokudo1[i,99],abs(seki_transfer1[i,100]),s=10,c='red', marker='.')
ax1.scatter(gunsokudo1[len(sekinums)-1,99],abs(seki_transfer1[len(sekinums)-1,100]),label='1.0MHz',s=10,c='red', marker='.')
for i in range(len(sekinums)-1):
	ax1.scatter(gunsokudo1[i,149],abs(seki_transfer1[i,150]),s=10,c='g', marker='.')
ax1.scatter(gunsokudo1[len(sekinums)-1,149],abs(seki_transfer1[len(sekinums)-1,150]),label='1.5MHz',s=10,c='g', marker='.')
#ax1.scatter(mean_gunsokudo1[49:149],abs(mean_transferx[50:150]),c='y',label='mean transfer', marker='.')
ax1.scatter(mean_gunsokudo4[49],sekit[50], s=200, c='y',label='mean_0.5MHz',marker='.')
ax1.scatter(mean_gunsokudo4[99],sekit[100], s=200, c=(0.3, 0.2, 0.8) ,label='mean_1.0MHz',marker='.')
ax1.scatter(mean_gunsokudo4[149],sekit[150],s=200, c='m', label='mean_1.5MHz',marker='.')
ax1.set_title('Covariance')
ax1.set_xlabel('group_vel[km/s]')
ax1.set_ylabel('amplitude')
ax1.set_xlim([0,10])
ax1.set_ylim([0,1.0])
ax1.grid(True)
ax1.legend()
fig1.savefig("kyoubunsanvel_seki.png",bbox_inches="tight")

fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
for i in range(len(karinums)-1):
	ax2.scatter(gunsokudo2[i,49],abs(kari_transfer1[i,50]),s=10, c='b',marker='.')
ax2.scatter(gunsokudo2[len(karinums)-1,49],abs(kari_transfer1[len(karinums)-1,50]),label='0.5MHz', s=10,c='b',marker='.')
for i in range(len(karinums)-1):
	ax2.scatter(gunsokudo2[i,99],abs(kari_transfer1[i,100]),s=10,c='red', marker='.')
ax2.scatter(gunsokudo2[len(karinums)-1,99],abs(kari_transfer1[len(karinums)-1,100]),label='1.0MHz',s=10,c='red', marker='.')
for i in range(len(karinums)-1):
	ax2.scatter(gunsokudo2[i,149],abs(kari_transfer1[i,150]),s=10,c='g', marker='.')
ax2.scatter(gunsokudo2[len(karinums)-1,149],abs(kari_transfer1[len(karinums)-1,150]),label='1.5MHz',s=10,c='g', marker='.')
#ax2.scatter(mean_gunsokudo2[49:149],abs(mean_transfery[50:150]),c='y', label='mean transfer',marker='.')
ax2.scatter(mean_gunsokudo5[49],karit[50],s=200,c='y',label='mean_0.5MHz', marker='.')
ax2.scatter(mean_gunsokudo5[99],karit[100],s=200,c=(0.3, 0.2, 0.8) ,label='mean_1.0MHz', marker='.')
ax2.scatter(mean_gunsokudo5[149],karit[150],s=200,c='m',label='mean_1.5MHz', marker='.')
ax2.set_title('Covariance')
ax2.set_xlabel('group_vel[km/s]')
ax2.set_ylabel('amplitude')
ax2.set_xlim([0,10])
ax2.set_ylim([0,1.0])
ax2.grid(True)
ax2.legend()
fig2.savefig("kyoubunsanvel_kari.png",bbox_inches="tight")


fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
for i in range(len(natonums)-1):
	ax3.scatter(gunsokudo3[i,49],abs(nato_transfer1[i,50]), s=10,c='b',marker='.')
ax3.scatter(gunsokudo3[len(natonums)-1,49],abs(nato_transfer1[len(natonums)-1,50]),label='0.5MHz', s=10,c='b',marker='.')
for i in range(len(natonums)-1):
	ax3.scatter(gunsokudo3[i,99],abs(nato_transfer1[i,100]),s=10,c='red', marker='.')
ax3.scatter(gunsokudo3[len(natonums)-1,99],abs(nato_transfer1[len(natonums)-1,100]),label='1.0MHz',s=10,c='red', marker='.')
for i in range(len(natonums)-1):
	ax3.scatter(gunsokudo3[i,149],abs(nato_transfer1[i,150]),s=10,c='g', marker='.')
ax3.scatter(gunsokudo3[len(natonums)-1,149],abs(nato_transfer1[len(natonums)-1,150]),label='1.5MHz',s=10,c='g', marker='.')
#ax3.scatter(mean_gunsokudo3[49:149],abs(mean_transferz[50:150]),c='y',label='mean transfer', marker='.')
ax3.scatter(mean_gunsokudo6[49],natot[50],s=200,c='y',label='mean_0.5MHz', marker='.')
ax3.scatter(mean_gunsokudo6[99],natot[100],s=200,c=(0.3, 0.2, 0.8) ,label='mean_1.0MHz', marker='.')
ax3.scatter(mean_gunsokudo6[149],natot[150],s=200,c='m',label='mean_1.5MHz', marker='.')
ax3.set_title('Covariance')
ax3.set_xlabel('group_vel[km/s]')
ax3.set_ylabel('amplitude')
ax3.set_xlim([0,10])
ax3.set_ylim([0,1.0])
ax3.grid(True)
ax3.legend()
fig3.savefig("kyoubunsanvel_nato.png",bbox_inches="tight")





#------------  散布図 振幅　群遅延直線のタイプ  -----------
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
for i in range(len(sekinums)):
	ax1.scatter(guntien1[i,149],abs(seki_transfer1[i,150]),c='b', marker='.')
ax1.scatter(mean_guntien1[149],abs(mean_transferx[150]),c='red', marker='.')
ax1.set_title('Covariance')
ax1.set_xlabel('group_delay[μsec]')
ax1.set_ylabel('amplitude')
ax1.set_xlim([0,2])
ax1.set_ylim([0,1.3])
ax1.grid(True)
fig1.savefig("kyoubunsan_150seki.png",bbox_inches="tight")

fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
for i in range(len(karinums)):
	ax2.scatter(guntien2[i,149],abs(kari_transfer1[i,150]),c='b', marker='.')
ax2.scatter(mean_guntien2[149],abs(mean_transfery[150]),c='red', marker='.')
ax2.set_title('Covariance')
ax2.set_xlabel('group_delay[μsec]')
ax2.set_ylabel('amplitude')
ax2.set_xlim([0,2])
ax2.set_ylim([0,1.3])
ax2.grid(True)
fig2.savefig("kyoubunsan_150kari.png",bbox_inches="tight")

fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
for i in range(len(natonums)):
	ax3.scatter(guntien3[i,149],abs(nato_transfer1[i,150]),c='b', marker='.')
ax3.scatter(mean_guntien3[149],abs(mean_transferz[150]),c='r', marker='.')
ax3.set_title('Covariance')
ax3.set_xlabel('group_delay[μsec]')
ax3.set_ylabel('amplitude')
ax3.set_xlim([0,2])
ax3.set_ylim([0,1.3])
ax3.grid(True)
fig3.savefig("kyoubunsan_150nato.png",bbox_inches="tight")





#------------  群速度 show python program  -----------
fig1=plt.figure()
ax=fig1.add_subplot(111)
for i in range(len(sekinums)):
	ax.plot(freq1,gunsokudo1[i,:],'b',label="<$trans$>")
ax.plot(freq1,mean_gunsokudo1,'r',label="<$trans$>")
ax.grid(True)

fsz=10; # fontsize
ax.set_xlabel("frequency(MHz)",fontsize=fsz)
ax.set_ylabel("group_vel[km/s]",fontsize=fsz)

ax.tick_params(labelsize=fsz)
ax.set_xlim([0.3,1.7])
ax.set_ylim([0,15])
fig1.savefig("gunsokudo_allseki.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
for i in range(len(karinums)):
	ax.plot(freq1,gunsokudo2[i,:],'b',label="<$trans$>")
ax.plot(freq1,mean_gunsokudo2,'r',label="<$trans$>")
ax.grid(True)

fsz=10; # fontsize
ax.set_xlabel("frequency(MHz)",fontsize=fsz)
ax.set_ylabel("group_vel[km/s]",fontsize=fsz)

ax.tick_params(labelsize=fsz)
ax.set_xlim([0.3,1.7])
ax.set_ylim([0,15])
fig1.savefig("gunsokudo_allkari.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
for i in range(len(natonums)):
	ax.plot(freq1,gunsokudo3[i,:],'b',label="<$trans$>")
ax.plot(freq1,mean_gunsokudo3,'r',label="<$trans$>")
ax.grid(True)

fsz=10; # fontsize
ax.set_xlabel("frequency(MHz)",fontsize=fsz)
ax.set_ylabel("group_vel[km/s]",fontsize=fsz)

ax.tick_params(labelsize=fsz)
ax.set_xlim([0.3,1.7])
ax.set_ylim([0,15])
fig1.savefig("gunsokudo_mean_nato.png",bbox_inches="tight")

#------------  群遅延 show python program  -----------
fig1=plt.figure()
ax=fig1.add_subplot(111)
for i in range(len(sekinums)):
	ax.plot(freq1,guntien1[i,:],'b',label="<$trans$>")
ax.plot(freq1,mean_guntien1,'r',label="<$trans$>")
ax.grid(True)

fsz=10; # fontsize
ax.set_xlabel("frequency(MHz)",fontsize=fsz)
ax.set_ylabel("group_delay[μsec]",fontsize=fsz)

ax.tick_params(labelsize=fsz)
ax.set_xlim([0.3,1.7])
ax.set_ylim([0,2])
fig1.savefig("guntien_allseki.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
for i in range(len(karinums)):
	ax.plot(freq1,guntien2[i,:],'b',label="<$trans$>")
ax.plot(freq1,mean_guntien2,'r',label="<$trans$>")
ax.grid(True)

fsz=10; # fontsize
ax.set_xlabel("frequency(MHz)",fontsize=fsz)
ax.set_ylabel("group_delay[μsec]",fontsize=fsz)

ax.tick_params(labelsize=fsz)
ax.set_xlim([0.3,1.7])
ax.set_ylim([0,2])
fig1.savefig("guntien_allkari.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
for i in range(len(natonums)):
	ax.plot(freq1,guntien3[i,:],'b',label="<$trans$>")
ax.plot(freq1,mean_guntien3,'r',label="<$trans$>")
ax.grid(True)

fsz=10; # fontsize
ax.set_xlabel("frequency(MHz)",fontsize=fsz)
ax.set_ylabel("group_delay[μsec]",fontsize=fsz)

ax.tick_params(labelsize=fsz)
ax.set_xlim([0.3,1.7])
ax.set_ylim([0,2])
fig1.savefig("guntien_allnato.png",bbox_inches="tight")


for i in range(Nx):
	fig1=plt.figure()
	ax=fig1.add_subplot(111)
	#ax.set_yscale('log')
	ax.plot(freq[:nf2],np.abs(declination[i,:nf2]),label="<$trans$>")
	ax.plot(freq[:nf2],np.abs(mean_declination[:nf2]),label="<$trans$>")
	ax.grid(True)

	fsz=10; # fontsize
	ax.set_xlabel("frequency(KHz)",fontsize=fsz)
	ax.set_ylabel("phase angle[rad]",fontsize=fsz)

	ax.tick_params(labelsize=fsz)
	#ax.set_xlim([0,10000])
	#ax.set_ylim([0,70])
	fig1.savefig("declination"+str(i)+".png",bbox_inches="tight")

	#ax.legend()


fig1=plt.figure()
ax=fig1.add_subplot(111)
ax.plot(time,mean_amp_original)
ax.plot(time,mean_amp)
ax.grid(True)
fsz=10; # fontsize
#ax.set_title("amplitude" ,fontsize=fsz)
ax.set_xlabel("time  (μsec)",fontsize=fsz)
ax.set_ylabel("amplitude",fontsize=fsz)

#ax.tick_params(labelsize=fsz)
ax.set_xlim([time[0],time[-1]])
ax.set_ylim([-0.04,0.04])
fig1.savefig("mean_amp.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
ax.plot(freq,abs(mean_Amp_original))
ax.plot(freq,abs(mean_Amp))
ax.grid(True)
fsz=10; # fontsize
#ax.set_title("amplitude" ,fontsize=fsz)
ax.set_xlabel("frequency (MHz)",fontsize=fsz)
ax.set_ylabel("Fourier amplitude",fontsize=fsz)
#ax.set_yscale('log')
#ax.tick_params(labelsize=fsz)
ax.set_xlim([0,2.5])
ax.set_ylim([0,1])
#ax.set_ylim([1.e-04,1]);
fig1.savefig("mean_Amp.png",bbox_inches="tight")



for i in range(Nx):
	fig1=plt.figure()
	ax=fig1.add_subplot(111)
	ax.plot(time,amp_re[i])
	ax.grid(True)
	fsz=10; # fontsize
	ax.set_title("Transfer function" + str(i),fontsize=fsz)
	ax.set_xlabel("time  (μsec)",fontsize=fsz)
	ax.set_ylabel("Transfer function",fontsize=fsz)

	ax.tick_params(labelsize=fsz)
	#ax.set_xlim([0,2])
	ax.set_ylim([-4.0,4.0])
	fig1.savefig("time_trans"+str(i)+".png",bbox_inches="tight")

	#ax.legend()

#------------  transfer show python program  -----------
fig1=plt.figure()
ax=fig1.add_subplot(111)
for i in range(len(sekinums)-1):
	ax.plot(freq,abs(seki_transfer1[i,:]),'b')
ax.plot(freq,abs(seki_transfer1[len(sekinums)-1,:]),'b',label="transfer function")
ax.plot(freq,np.std(abs(seki_transfer1),0),'m',label="std")
ax.plot(freq,abs(mean_transferx),'r',label="mean")
ax.grid(True)
ax.set_xlim([0,2.0])
ax.set_ylim([0,1.0])
#ax.set_ylim([0,1])
ax.set_title("seki_transfer")
ax.set_xlabel("frequency (MHz)");
ax.set_ylabel("transfer function");
ax.legend()
fig1.savefig("seki_transfer.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
for i in range(len(karinums)-1):
	ax.plot(freq,abs(kari_transfer1[i,:]),'b')
ax.plot(freq,abs(kari_transfer1[len(karinums)-1,:]),'b',label="transfer function")
ax.plot(freq,np.std(abs(kari_transfer1),0),'m',label="std")
ax.plot(freq,abs(mean_transfery),'r',label="mean")
ax.grid(True)
ax.set_xlim([0,2.0])
ax.set_ylim([0,1.0])
#ax.set_ylim([0,1])
ax.set_title("kari_transfer")
ax.set_xlabel("frequency (MHz)");
ax.set_ylabel("transfer function");
ax.legend()
fig1.savefig("kari_transfer.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
for i in range(len(natonums)-1):
	ax.plot(freq,abs(nato_transfer1[i,:]),'b')
ax.plot(freq,abs(nato_transfer1[len(natonums)-1,:]),'b',label="transfer function")
ax.plot(freq,np.std(abs(nato_transfer1),0),'m',label="std")
ax.plot(freq,abs(mean_transferz),'r',label="mean")
ax.grid(True)
ax.set_xlim([0,2.0])
ax.set_ylim([0,1.0])
#ax.set_ylim([0,1])
ax.set_title("nato_transfer")
ax.set_xlabel("frequency (MHz)");
ax.set_ylabel("transfer function");
ax.legend()
fig1.savefig("nato_transfer.png",bbox_inches="tight")



#------------  位相 show python program  -----------
fig1=plt.figure()
ax=fig1.add_subplot(111)
for i in range(len(sekinums)):
	ax.plot(freq,seki_declination[i,:],'b')
#ax.plot(freq,-np.std(seki_declination,0),'b')
ax.plot(freq,seki_mean_declination,'r')
ax.grid(True)
ax.set_xlim([0,2.0])
#ax.set_ylim([5,-20])
ax.set_ylim([2,-20])
ax.set_title("mean_seki_declination")
ax.set_xlabel("frequency (MHz)");
ax.set_ylabel("phase angle [rad]");
fig1.savefig("seki_std_declination12.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
for i in range(len(karinums)):
	ax.plot(freq,kari_declination[i,:],'b')
#ax.plot(freq,-np.std(kari_declination,0),'b')
ax.plot(freq,kari_mean_declination,'r')
ax.grid(True)
ax.set_xlim([0,2.0])
#ax.set_ylim([5,-20])
ax.set_ylim([2,-20])
ax.set_title("mean_kari_declination")
ax.set_xlabel("frequency (MHz)");
ax.set_ylabel("phase angle [rad]");
fig1.savefig("kari_std_declination12.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
for i in range(len(natonums)):
	ax.plot(freq,nato_declination[i,:],'b')
#ax.plot(freq,-np.std(nato_declination,0),'b')
ax.plot(freq,nato_mean_declination,'r')
ax.grid(True)
ax.set_xlim([0,2.0])
#ax.set_ylim([5,-20])
ax.set_ylim([2,-20])
ax.set_title("mean_nato_declination")
ax.set_xlabel("frequency (MHz)");
ax.set_ylabel("phase angle [rad]");
fig1.savefig("nato_std_declination12.png",bbox_inches="tight")

'''
#------------  群遅延 show python program  -----------

fig1=plt.figure()
ax=fig1.add_subplot(111)
results, edges = np.histogram(flatguntien_seki, bins=100, range=(0,2), normed=True)
binWidth = edges[1] - edges[0]
ax.bar(edges[:-1], results*binWidth, binWidth)
#ax.hist(seki_guntien,bins=100, density=True)
ax.set_xlim([0,2])
ax.set_ylim([0,0.06])
ax.grid(True)
fig1.savefig("seki_tienall_hist.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
results, edges = np.histogram(flatguntien_kari, bins=100, range=(0,2), normed=True)
binWidth = edges[1] - edges[0]
ax.bar(edges[:-1], results*binWidth, binWidth)
#ax.hist(kari_guntien,bins=100, density=True)
ax.set_xlim([0,2])
ax.set_ylim([0,0.06])
ax.grid(True)
fig1.savefig("kari_tienall_hist.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
results, edges = np.histogram(flatguntien_nato, bins=100, range=(0,2), normed=True)
binWidth = edges[1] - edges[0]
ax.bar(edges[:-1], results*binWidth, binWidth)
#ax.hist(seki_guntien,bins=100, density=True)
ax.set_xlim([0,2])
ax.set_ylim([0,0.06])
ax.grid(True)
fig1.savefig("nato_tienall_hist.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
results, edges = np.histogram(flatguntien_seki, bins=100, range=(0,2), normed=True)
binWidth = edges[1] - edges[0]
ax.bar(edges[:-1], results*binWidth, binWidth,label="quartz",alpha=0.5)
results, edges = np.histogram(flatguntien_kari, bins=100, range=(0,2), normed=True)
binWidth = edges[1] - edges[0]
ax.bar(edges[:-1], results*binWidth, binWidth, label="potassium feldspar",alpha=0.5)
results, edges = np.histogram(flatguntien_nato, bins=100, range=(0,2), normed=True)
binWidth = edges[1] - edges[0]
ax.bar(edges[:-1], results*binWidth, binWidth, label="plagioclase",alpha=0.5)
ax.set_xlim([0,2])
ax.set_ylim([0,0.06])
ax.set_xlabel("group delay[μsec]");
ax.grid(True)
ax.legend()
fig1.savefig("tien_hist.png",bbox_inches="tight")

'''

#------------  群速度 show python program  -----------

fig1=plt.figure()
ax=fig1.add_subplot(111)
results, edges = np.histogram(flatgunsokudo_seki, bins=50, range=(0,10), normed=True)
binWidth = edges[1] - edges[0]
ax.bar(edges[:-1], results*binWidth, binWidth)
#ax.hist(seki_guntien,bins=100, density=True)
#ax.set_xlim([0,10])
ax.set_ylim([0,0.1])
ax.grid(True)
fig1.savefig("seki_gunsokudoall_hist.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
results, edges = np.histogram(flatgunsokudo_kari, bins=50, range=(0,10), normed=True)
binWidth = edges[1] - edges[0]
ax.bar(edges[:-1], results*binWidth, binWidth)
#ax.hist(kari_guntien,bins=100, density=True)
#ax.set_xlim([0,10])
ax.set_ylim([0,0.1])
ax.grid(True)
fig1.savefig("kari_gunsokudoall_hist.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
results, edges = np.histogram(flatgunsokudo_nato, bins=50, range=(0,10), normed=True)
binWidth = edges[1] - edges[0]
ax.bar(edges[:-1], results*binWidth, binWidth)
#ax.hist(seki_guntien,bins=100, density=True)
#ax.set_xlim([0,10])
ax.set_ylim([0,0.1])
ax.grid(True)
fig1.savefig("nato_gunsokudoall_hist.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
results, edges = np.histogram(flatgunsokudo_seki, bins=50, range=(0,10), normed=True)
binWidth = edges[1] - edges[0]
ax.bar(edges[:-1], results*binWidth, binWidth,label="seki",alpha=0.5)
results, edges = np.histogram(flatgunsokudo_kari, bins=50, range=(0,10), normed=True)
binWidth = edges[1] - edges[0]
ax.bar(edges[:-1], results*binWidth, binWidth, label="kari",alpha=0.5)
results, edges = np.histogram(flatgunsokudo_nato, bins=50, range=(0,10), normed=True)
binWidth = edges[1] - edges[0]
ax.bar(edges[:-1], results*binWidth, binWidth, label="nato",alpha=0.5)
#ax.set_xlim([0,10])
ax.set_ylim([0,0.1])
ax.set_xlabel("group velocity[km/sec]");
ax.grid(True)
ax.legend()
fig1.savefig("gunsokudo_hist.png",bbox_inches="tight")



#------------  transfer probability show  -----------------------

fig1=plt.figure()
ax=fig1.add_subplot(111, projection="3d") # 3Dの軸を作成
relative = np.ones_like(z_seki[0,:])/len(z_seki[0,:])
hist, xedges, yedges = np.histogram2d(z_seki[0,:], z_seki[1,:], bins=100, weights=relative) # 100個区切りでxとyのヒストグラムを作成
xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1]) # x,y座標を3D用の形式に変換（その１）
zpos = 0 # zは常に0を始点にする
dx = abs(xpos[0][1]) - abs(xpos[0][0]) # x座標の幅を設定
dy = abs(ypos[1][0]) - abs(ypos[0][0]) # y座標の幅を設定
dz = hist.ravel() # z座標の幅は棒の長さに相当

xpos = xpos.ravel() # x座標を3D用の形式に変換（その２）
ypos = ypos.ravel() # y座標を3D用の形式に変換（その２）

ax.bar3d(xpos,ypos,zpos,dx,dy,dz) # ヒストグラムを3D空間に表示
plt.title("Histogram 2D") # タイトル表示
plt.xlabel("tilt") # x軸の内容表示
plt.ylabel("max") # y軸の内容表示
ax.set_zlabel("Z") # z軸の内容表示
ax.set_xlim([-0.5,0.1])
ax.set_ylim([0.1,0.6])
ax.set_zlim([0,0.1])
ax.grid(True)
fig1.savefig("seki_transfer_probability.png")

fig1=plt.figure()
ax=fig1.add_subplot(111, projection="3d") # 3Dの軸を作成
relative = np.ones_like(z_kari[0,:])/len(z_kari[0,:])
hist, xedges, yedges = np.histogram2d(z_kari[0,:], z_kari[1,:], bins=100, weights=relative) # 100個区切りでxとyのヒストグラムを作成
xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1]) # x,y座標を3D用の形式に変換（その１）
zpos = 0 # zは常に0を始点にする
dx = abs(xpos[0][1]) - abs(xpos[0][0]) # x座標の幅を設定
dy = abs(ypos[1][0]) - abs(ypos[0][0]) # y座標の幅を設定
dz = hist.ravel() # z座標の幅は棒の長さに相当

xpos = xpos.ravel() # x座標を3D用の形式に変換（その２）
ypos = ypos.ravel() # y座標を3D用の形式に変換（その２）

ax.bar3d(xpos,ypos,zpos,dx,dy,dz) # ヒストグラムを3D空間に表示
plt.title("Histogram 2D") # タイトル表示
plt.xlabel("tilt") # x軸の内容表示
plt.ylabel("max") # y軸の内容表示
ax.set_zlabel("Z") # z軸の内容表示
ax.set_xlim([-0.5,0.1])
ax.set_ylim([0.1,0.6])
ax.set_zlim([0,0.1])
ax.grid(True)
fig1.savefig("kari_transfer_probability.png")

fig1=plt.figure()
ax=fig1.add_subplot(111, projection="3d") # 3Dの軸を作成
relative = np.ones_like(z_nato[0,:])/len(z_nato[0,:])
hist, xedges, yedges = np.histogram2d(z_nato[0,:], z_nato[1,:], bins=100, weights=relative) # 100個区切りでxとyのヒストグラムを作成
xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1]) # x,y座標を3D用の形式に変換（その１）
zpos = 0 # zは常に0を始点にする
dx = abs(xpos[0][1]) - abs(xpos[0][0]) # x座標の幅を設定
dy = abs(ypos[1][0]) - abs(ypos[0][0]) # y座標の幅を設定
dz = hist.ravel() # z座標の幅は棒の長さに相当

xpos = xpos.ravel() # x座標を3D用の形式に変換（その２）
ypos = ypos.ravel() # y座標を3D用の形式に変換（その２）

ax.bar3d(xpos,ypos,zpos,dx,dy,dz) # ヒストグラムを3D空間に表示
plt.title("Histogram 2D") # タイトル表示
plt.xlabel("tilt") # x軸の内容表示
plt.ylabel("max") # y軸の内容表示
ax.set_zlabel("Z") # z軸の内容表示
ax.set_xlim([-0.5,0.1])
ax.set_ylim([0.1,0.6])
ax.set_zlim([0,0.1])
ax.grid(True)
fig1.savefig("nato_transfer_probability.png")



fig1=plt.figure()
ax=fig1.add_subplot(111)
results, edges = np.histogram(seki_transfermax, bins=100, normed=True)
binWidth = abs(edges[1] - edges[0])
ax.bar(edges[:-1], results*binWidth, binWidth)
#ax.hist(kari_guntien,bins=100, density=True)
#ax.set_xlim([0,10])
ax.set_ylim([0,0.1])
ax.grid(True)
fig1.savefig("kari_transfer_max.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
results, edges = np.histogram(kari_transfermax, bins=100, normed=True)
binWidth = abs(edges[1] - edges[0])
ax.bar(edges[:-1], results*binWidth, binWidth)
#ax.hist(kari_guntien,bins=100, density=True)
#ax.set_xlim([0,10])
ax.set_ylim([0,0.1])
ax.grid(True)
fig1.savefig("kari_transfermax.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
results, edges = np.histogram(nato_transfermax, bins=100, normed=True)
binWidth = abs(edges[1] - edges[0])
ax.bar(edges[:-1], results*binWidth, binWidth)
#ax.hist(seki_guntien,bins=100, density=True)
#ax.set_xlim([0,10])
ax.set_ylim([0,0.1])
ax.grid(True)
fig1.savefig("nato_transfermax.png",bbox_inches="tight")

fig1=plt.figure()
ax=fig1.add_subplot(111)
results, edges = np.histogram(seki_transfermax, bins=100, normed=True)
binWidth = abs(edges[1] - edges[0])
ax.bar(edges[:-1], results*binWidth, binWidth,label="seki",alpha=0.5)
results, edges = np.histogram(kari_transfermax, bins=100, normed=True)
binWidth = abs(edges[1] - edges[0])
ax.bar(edges[:-1], results*binWidth, binWidth, label="kari",alpha=0.5)
results, edges = np.histogram(nato_transfermax, bins=100, normed=True)
binWidth = abs(edges[1] - edges[0])
ax.bar(edges[:-1], results*binWidth, binWidth, label="nato",alpha=0.5)
#ax.set_xlim([0,10])
ax.set_ylim([0,0.1])
ax.grid(True)
ax.legend()
fig1.savefig("transfer_max.png",bbox_inches="tight")

fig = plt.figure()
ax = fig.add_subplot(111)
relative1 = np.ones_like(z_seki[0,:])/len(z_seki[0,:])
relative2 = np.ones_like(z_seki[1,:])/len(z_seki[1,:])
H = ax.hist2d(z_seki[0,:],z_seki[1,:], bins=100, weights=(relative1,relative2), cmap=cm.jet)
ax.set_title('4th graph')
ax.set_xlabel('x')
ax.set_ylabel('y')

fig.colorbar(H[3],ax=ax)
plt.show()


#------------  transfer function 3-D show  -----------------------

fig1=plt.figure()
ax=fig1.add_subplot(111, projection="3d") # 3Dの軸を作成
relative = np.ones_like(z_nato[0,:])/len(z_nato[0,:])
hist, xedges, yedges = np.histogram2d(z_seki[0:len(sekinums),:], z_seki[len(sekinums),:], bins=100, weights=relative) # 100個区切りでxとyのヒストグラムを作成
xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1]) # x,y座標を3D用の形式に変換（その１）
zpos = 0 # zは常に0を始点にする
dx = abs(xpos[0][1]) - abs(xpos[0][0]) # x座標の幅を設定
dy = abs(ypos[1][0]) - abs(ypos[0][0]) # y座標の幅を設定
dz = hist.ravel() # z座標の幅は棒の長さに相当

xpos = xpos.ravel() # x座標を3D用の形式に変換（その２）
ypos = ypos.ravel() # y座標を3D用の形式に変換（その２）

ax.bar3d(xpos,ypos,zpos,dx,dy,dz) # ヒストグラムを3D空間に表示
plt.title("Histogram 2D") # タイトル表示
plt.xlabel("tilt") # x軸の内容表示
plt.ylabel("max") # y軸の内容表示
ax.set_zlabel("Z") # z軸の内容表示
ax.set_xlim([-0.5,0.1])
ax.set_ylim([0.1,0.6])
ax.set_zlim([0,0.1])
ax.grid(True)
fig1.savefig("seki_transfer_probability.png")

fig, ax = plt.subplots()
im=ax.imshow(seki_count,vmin=0,vmax=1,cmap="jet")
#heatmap = ax.pcolor(seki_count, cmap=plt.cm.Blues,vmax=1,vmin=0)
#ax.set_ylim(50,150)
plt.show()
'''
