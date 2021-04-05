import numpy as np
import matplotlib.pyplot as plt
import Ascan as Asc

if __name__=="__main__":

    dir_name="./"
    fname=dir_name+"/1MHznew.csv"

    dir_name="./K_Feldspar"
    num=329
    fname=dir_name+"/scope_"+str(num)+".csv"

    awv=Asc.AWV()
    awv.load(fname)

    fig1=plt.figure()
    ax=fig1.add_subplot(111)
    awv.plot_ascan(ax)


    tb=12.3; w_6dB=1;mexp=6
    tb=13.0; w_6dB=1;mexp=6

    awv.Butterworth(tb,w_6dB,6,apply=True)
    awv.plot_ascan(ax)

    fig2=plt.figure()
    bx1=fig2.add_subplot(211)
    bx2=fig2.add_subplot(212)
    awv.plot_FFT(bx1)
    awv.plot_gdelay(bx2)


    plt.show()
