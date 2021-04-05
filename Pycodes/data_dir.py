class DATA_DIR:
    def __init__(self,Mineral):
        self.Mineral=Mineral
        if Mineral=="Na":
            dir_name="Na_Feldspar"
            tb=12.5; tw_6dB=1; mexp=6
            nfile=548;
            bT=2.03 # Na
            rho=2.62
        elif Mineral=="K":
            dir_name="K_Feldspar"
            tb=12.5; tw_6dB=1; mexp=6
            nfile=824;
            bT=1.89 # K
            rho=2.56
        elif Mineral=="Qt":
            dir_name="Quartz"
            tb=12.5; tw_6dB=1.0; mexp=6
            nfile=589;
            bT=2.08 #Qt
            rho=2.65
        elif Mineral=="Ref":
            dir_name="./"
            tb=12; tw_6dB=1; mexp=6
            nfile=1;
            rho=1.05
            self.cp=2.4
            bT=1 
        else:
            print("no data for ",Mineral)
            exit();
        self.dir_name=dir_name
        self.tb=tb
        self.tw_6dB=tw_6dB
        self.mexp=mexp
        self.nfile=nfile
        self.rho=rho
        self.bT=bT
    def get_dir_name(self):
        return(self.dir_name)
    def get_win_params(self):
        return(self.tb, self.tw_6dB, self.mexp)
    def get_nfile(self):
        return(self.nfile)
    def get_bT(self):
        return(self.bT)
    def show(self):
        print("data directory= ",self.dir_name)
        print("win params(tb, tw, mexp)= ",self.tb,self.tw_6dB,self.mexp)
        print("number of files= ",self.nfile)
        print("efficiency factor= ",self.bT)

