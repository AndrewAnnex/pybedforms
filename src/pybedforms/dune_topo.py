import numpy as np
from numpy import pi, sin, cos
import numba as nb
from typing import Tuple


#jitclass
class DuneTopo(object):

    def __init__(self):
        self.SPCNGF = 50.0  # (3)  WAVELENGTH OF BEDFORMS IN FIRST SET
        self.PHASEF = 0.0  # (4)  Bedform phase (controls placement within block diagram)
        self.SMTRYF = 0.0  # (5)  Mean asymmetry
        self.SMCHGF = 0.0  # (6)  Amplitude of asymmetry cycle
        self.SMPRDF = 1.0  # (7)  Period of asymmetry cycle
        self.SMFAZF = 0.0  # (8)  Phase of asymmetry cycle
        self.HTRTOF = 1.0  # (9)  Mean steepness
        self.HTCHGF = 0.0  # (10) Amplitude of steepness cycle
        self.HTPRDF = 1.0  # (11) Period of steepness cycle
        self.HTFAZF = 0.0  # (12) Phase of steepness cycle
        self.SNSPF1 = 0.0  # (13) Wavelength (alon g -crest) of first set of pla n -form sinuosities
        self.SNMGF1 = 0.0  # (14) Amplitude (measured in pla n -form) of first set of sinuosities
        self.SNFZF1 = 0.0  # (15) Phase of first set of sinuosities
        self.SNVLF1 = 0.0  # (16) Migration speed (alon g -crest) of first set of sinuosities
        self.SNSPF2 = 0.0  # (17) Wavelength (alon g -crest) of second set of pla n -form sinuosities
        self.SNMGF2 = 0.0  # (18) Amplitude (measured in pla n -form) of second set of sinuosities
        self.SNFZF2 = 0.0  # (19) Phase of second set of sinuosities
        self.SNVLF2 = 0.0  # (20) Migration speed (alon g -crest) of second set of sinuosities
        self.TRENDF = 90.0  # (21) Migration direction of bedform
        self.VELOCF = 0.0  # (22) Mean migration speed of bedform
        self.VLCHGF = 0.0  # (23) Amplitude of speed cycle
        self.VLPRDF = 1.0  # (24) Period of speed cycle
        self.VLFAZF = 0.0  # (25) Phase of speed cycle
        self.SPCNGS = 0.0  # (26) WAVELENGTH OF BEDFORMS IN SECOND SET
        self.PHASES = 0.0  # (27) Bedform phase (controls placement within block diagram)
        self.SMTRYS = 0.0  # (28) Mean asymmetry
        self.SMCHGS = 0.0  # (29) Amplitude of asymmetry cycle
        self.SMPRDS = 1.0  # (30) Period of asymmetry cycle
        self.SMFAZS = 0.0  # (31) Phase of asymmetry cycle
        self.HTRTOS = 0.0  # (32) Mean steepness
        self.HTCHGS = 0.0  # (33) Amplitude of steepness cycle
        self.HTPRDS = 1.0  # (34) Period of steepness cycle
        self.HTFAZS = 0.0  # (35) Phase of steepness cycle
        self.SNSPS1 = 0.0  # (36) Wavelength (alon g -crest) of first set of pla n -form sinuosities
        self.SNMGS1 = 0.0  # (37) Amplitude (measured in pla n -form) of first set of sinuosities
        self.SNFZS1 = 0.0  # (38) Phase of first set of sinuosities
        self.SNVLS1 = 0.0  # (39) Migration speed (alon g -crest) of first set of sinuosities
        self.SNSPS2 = 0.0  # (40) Wavelength (alon g -crest) of second set of pla n -form sinuosities
        self.SNMGS2 = 0.0  # (41) Amplitude (measured in pla n -form) of second set of sinuosities
        self.SNFZS2 = 0.0  # (42) Phase of second set of sinuosities
        self.SNVLS2 = 0.0  # (43) Migration speed (alon g -crest) of second set of sinuosities
        self.TRENDS = 0.0  # (44) Migration direction of bedform
        self.VELOCS = 0.0  # (45) Mean migration speed of bedform
        self.VLCHGS = 1.0  # (46) Amplitude of speed cycle
        self.VLPRDS = 1.0  # (47) Period of speed cycle
        self.VLFAZS = 0.0  # (48) Phase of speed cycle
        self.SPCNGT = 0.0  # (49) WAVELENGTH OF BEDFORMS IN THIRD SET
        self.PHASET = 0.0  # (50) Bedform phase (controls placement within block diagram)
        self.SMTRYT = 0.0  # (51) Mean asymmetry
        self.SMCHGT = 0.0  # (52) Amplitude of asymmetry cycle
        self.SMPRDT = 1.0  # (53) Period of asymmetry cycle
        self.SMFAZT = 0.0  # (54) Phase of asymmetry cycle
        self.HTRTOT = 0.0  # (55) Mean steepness
        self.HTCHGT = 0.0  # (56) Amplitude of steepness cycle
        self.HTPRDT = 1.0  # (57) Period of steepness cycle
        self.HTFAZT = 0.0  # (58) Phase of steepness cycle
        self.SNSPT1 = 0.0  # (59) Wavelength (alon g -crest) of first set of pla n -form sinuosities
        self.SNMGT1 = 0.0  # (60) Amplitude (measured in pla n -form) of first set of sinuosities
        self.SNFZT1 = 0.0  # (61) Phase of first set of sinuosities
        self.SNVLT1 = 0.0  # (62) Migration speed (alon g -crest) of first set of sinuosities
        self.SNSPT2 = 0.0  # (63) Wavelength (alon g -crest) of second set of pla n -form sinuosities
        self.SNMGT2 = 0.0  # (64) Amplitude (measured in pla n -form) of second set of sinuosities
        self.SNFZT2 = 0.0  # (65) Phase of second set of sinuosities
        self.SNVLT2 = 0.0  # (66) Migration speed (alon g -crest) of second set of sinuosities
        self.TRENDT = 0.0  # (67) Migration direction of bedform
        self.VELOCT = 0.0  # (68) Mean migration speed of bedform
        self.VLCHGT = 0.0  # (69) Amplitude of speed cycle
        self.VLPRDT = 1.0  # (70) Period of speed cycle
        self.VLFAZT = 0.0  # (71) Phase of speed cycle
        self.TYPE = 2  # (72) Type of superpositioning (INTEGER)
        self.CHOICE = 0  # (73) Rotation option (INTEGER
        self.ELVMIN = -1.0  # (74) Elevation of interdune flats
        self.DEPRAT = .8  # (75) Rate of deposition
        self.DEPCHG = 0.0  # (76) Amplitude of cycle in rate of deposition
        self.DEPPRD = 1.0  # (77) Period of cycle in rate of deposition
        self.DEPFAZ = 0.0  # (78) Phase of cycle in rate of deposition
        self.NBEDSH = 2500  # (79) Time from t=0 to beginning of depositional episode (INTEGER)
        self.INTXBD = 4  # (80) Interval between drawing of crossbeds (INTEGER)
        self.FRMNUM = 1  # (81) Time from t=0 to end of depositional episode (INTEGER)
        self.CAPLOG = 'true'  # (82) Print caption? (LOGICAL VARIABLE)
        self.FILOG = 'true'  # (83) Print name of input file? (LOGICAL VARIABLE)
        self.FRMLOG = 'true'  # (84) Print time at end of depositional episode? (LOGICAL VARIABLE)
        self.IGRDSP = 1  # (85) Precision? (low number) or speed? (high number) (1 2 4 5 or 10)
        self.ZHORIZ = 1.0  # (86) Elevation of horizontal section


        # validation inits from dune init

        if self.SPCNGF == 0:
            self.FD = 1e5
        if self.SPCNGS == 0:
            self.SD = 1e5
        if self.SPCNGT == 0:
            self.TD = 1e5
        if self.SNSPF1 == 0:
            self.SNSPF1 = 1e5
        if self.SNSPF2 == 0:
            self.SNSPF2 = 1e5
        if self.SNSPS1 == 0:
            self.SNSPS1 = 1e5
        if self.SNSPS2 == 0:
            self.SNSPS2 = 1e5
        if self.SNSPT1 == 0:
            self.SNSPT1 = 1e5
        if self.SNSPT2 == 0:
            self.SNSPT2 = 1e5
        if self.SMPRDF == 0:
            self.SMPRDF = 1e5
        if self.SMPRDS == 0:
            self.SMPRDS = 1e5
        if self.SMPRDT == 0:
            self.SMPRDT = 1e5
        if self.HTPRDF == 0:
            self.HTPRDF = 1e5
        if self.HTPRDS == 0:
            self.HTPRDS = 1e5
        if self.HTPRDT == 0:
            self.HTPRDT = 1e5
        if self.VLPRDF == 0:
            self.VLPRDF = 1e5
        if self.VLPRDS == 0:
            self.VLPRDS = 1e5
        if self.VLPRDT == 0:
            self.VLPRDT = 1e5
        if self.DEPPRD == 0:
            self.DEPPRD = 1e5
        if self.SPCNGF != 0:
            self.FD = 100 / self.SPCNGF
        if self.SPCNGS != 0:
            self.SD = 100 / self.SPCNGS
        if self.SPCNGT != 0:
            self.TD = 100 / self.SPCNGT
        # Make sure ZHORIZ is not equal to 1.0.
        if np.abs(self.ZHORIZ - 1.0) < 0.001:
            self.ZHORIZ = 1.001


    def run(self, x, y, TIME, ZBED, run, GridSize, EdgeGridSpace):
        """
        main run method, below will be split into other functions as needed
        :return:
        """
        # Calculate deposition.
        DPOSIT = (self.DEPRAT*TIME)+self.DEPCHG*(-self.DEPPRD/(2*pi))*cos(((self.DEPFAZ-90)*2*pi/360)+(TIME*2*pi/self.DEPPRD))

        # Calculate bedform asymmetry.
        PROFF = ((1-self.SMTRYF)+self.SMCHGF*sin(((-TIME/self.SMPRDF)-self.SMFAZF/360)*pi*2.0))*pi/2.0
        PROFS = ((1-self.SMTRYS)+self.SMCHGS*sin(((-TIME/self.SMPRDS)-self.SMFAZS/360)*pi*2.0))*pi/2.0
        PROFT = ((1-self.SMTRYT)+self.SMCHGT*sin(((-TIME/self.SMPRDT)-self.SMFAZT/360)*pi*2.0))*pi/2.0

        # Calculate bedform sizes.
        sizef = self.HTRTOF+(self.HTCHGF*sin(((TIME/self.HTPRDF)*2*pi)+((self.HTFAZF-90)*2*pi/360.0)))
        sizes = self.HTRTOS+(self.HTCHGS*sin(((TIME/self.HTPRDS)*2*pi)+((self.HTFAZS-90)*2*pi/360.0)))
        sizet = self.HTRTOT+(self.HTCHGT*sin(((TIME/self.HTPRDT)*2*pi)+((self.HTFAZT-90)*2*pi/360.0)))

        # Calculate bedform locations.
        dispfd = self.VELOCF*TIME+self.VLCHGF*(-self.VLPRDF/(2*pi))*cos(((self.VLFAZF-90)*2*pi/360)+(TIME*2*pi/self.VLPRDF))
        dispsd = self.VELOCS*TIME+self.VLCHGS*(-self.VLPRDS/(2*pi))*cos(((self.VLFAZS-90)*2*pi/360)+(TIME*2*pi/self.VLPRDS))
        disptd = self.VELOCT*TIME+self.VLCHGT*(-self.VLPRDT/(2*pi))*cos(((self.VLFAZT-90)*2*pi/360)+(TIME*2*pi/self.VLPRDT))

        yfd = x*sin(self.TRENDF*2*pi/360)+y*cos(self.TRENDF*2*pi/360)
        ysd = x*sin(self.TRENDS*2*pi/360)+y*cos(self.TRENDS*2*pi/360)
        ytd = x*sin(self.TRENDT*2*pi/360)+y*cos(self.TRENDT*2*pi/360)
        xfd = (x*cos(self.TRENDF*2*pi/360)-y*sin(self.TRENDF*2*pi/360) - dispfd-self.SPCNGF*self.PHASEF/360) - \
              (self.SNMGF1*sin((yfd*2*pi/self.SNSPF1) + ((self.SNFZF1*2*pi/360) + (TIME*self.SNVLF1*2*pi/self.SNSPF1)))) - \
              (self.SNMGF2*sin((yfd*2*pi/self.SNSPF2) + ((self.SNFZF2*2*pi/360) + (TIME*self.SNVLF2*2*pi/self.SNSPF2))))
        xsd = (x*cos(self.TRENDS*2*pi/360)-y*sin(self.TRENDS*2*pi/360) - dispsd-self.SPCNGS*self.PHASES/360) - \
              (self.SNMGS1*sin((ysd*2*pi/self.SNSPS1) + ((self.SNFZS1*2*pi/360) + (TIME*self.SNVLS1*2*pi/self.SNSPS1)))) - \
              (self.SNMGS2*sin((ysd*2*pi/self.SNSPS2) + ((self.SNFZS2*2*pi/360) + (TIME*self.SNVLS2*2*pi/self.SNSPS2))))
        xtd = (x*cos(self.TRENDT*2*pi/360)-y*sin(self.TRENDT*2*pi/360) - disptd-self.SPCNGT*self.PHASET/360) - \
              (self.SNMGT1*sin((ytd*2*pi/self.SNSPT1) + ((self.SNFZT1*2*pi/360) + (TIME*self.SNVLT1*2*pi/self.SNSPT1)))) - \
              (self.SNMGT2*sin((ytd*2*pi/self.SNSPT2) + ((self.SNFZT2*2*pi/360) + (TIME*self.SNVLT2*2*pi/self.SNSPT2))))
        # Superimpose sets of dunes
        FD = self.FD
        SD = self.SD
        TD = self.TD
        z = None
        if self.TYPE == 1:
            zfd = (-6*sin(xfd*FD*pi/50)/FD-1.5*sin((xfd*FD*pi/25)+PROFF)/FD)
            shape = 1
            zfd = zfd*sizef
            zsd = (-6*sin(xsd*SD*pi/50)/SD-1.5*sin((xsd*SD*pi/25)+PROFS)/SD)*shape*sizes
            ztd = (-6*sin(xtd*TD*pi/50)/TD-1.5*sin((xtd*TD*pi/25)+PROFT)/TD)*shape*sizet
            z = zfd + zsd + ztd + DPOSIT
            z = max(z, ((7.5/FD)+(7.5/SD)+(7.5/TD))*self.ELVMIN+DPOSIT)
        elif self.TYPE == 2:
            zfd = (-6*sin(xfd*FD*pi/50)/FD-1.5*sin((xfd*FD*pi/25)+PROFF)/FD)
            shape = ((7.5/FD)-zfd)/(15/FD)
            zfd = zfd*sizef
            zsd = (-6*sin(xsd*SD*pi/50)/SD-1.5*sin((xsd*SD*pi/25)+PROFS)/SD)*shape*sizes
            ztd = (-6*sin(xtd*TD*pi/50)/TD-1.5*sin((xtd*TD*pi/25)+PROFT)/TD)*shape*sizet
            z = zfd + zsd + ztd + DPOSIT
            z = max (z, ((7.5/FD)+(7.5/SD)+(7.5/TD))*self.ELVMIN+DPOSIT)
        elif self.TYPE == 3:
            BD = 100 / max(self.SPCNGF, self.SPCNGS, self.SPCNGT)
            zfd = (-6*sin(xfd*FD*pi/50)/FD-1.5*sin((xfd*FD*pi/25)+PROFF)/FD+7.5/FD)*sizef
            zsd = (-6*sin(xsd*SD*pi/50)/SD-1.5*sin((xsd*SD*pi/25)+PROFS)/SD+7.5/SD)*sizes
            ztd = (-6*sin(xtd*TD*pi/50)/TD-1.5*sin((xtd*TD*pi/25)+PROFT)/TD+7.5/TD)*sizet
            z = max(zfd, zsd)
            z = max(z, ztd)
            z = z + DPOSIT
            z = max(z, (7.5*(1+self.ELVMIN) / BD + DPOSIT) )
        elif self.TYPE == 4:
            zfd = (-6*sin(xfd*FD*pi/50)/FD-1.5*sin((xfd*FD*pi/25)+PROFF)/FD)
            shape = 1-((7.5/FD)-zfd)/(15/FD)
            zfd = zfd*sizef
            zsd = (-6*sin(xsd*SD*pi/50)/SD-1.5*sin((xsd*SD*pi/25)+PROFS)/SD)*shape*sizes
            ztd = (-6*sin(xtd*TD*pi/50)/TD-1.5*sin((xtd*TD*pi/25)+PROFT)/TD)*shape*sizet
            z = zfd + zsd + ztd + DPOSIT
            z = max(z, ((7.5/FD)+(7.5/SD)+(7.5/TD))*self.ELVMIN+DPOSIT)
        elif self.TYPE == 5:
            zfd = (-6*sin(xfd*FD*pi/50)/FD-1.5*sin((xfd*FD*pi/25)+PROFF)/FD)*sizef
            zsd = (-6*sin(xsd*SD*pi/50)/SD-1.5*sin((xsd*SD*pi/25)+PROFS)/SD+7.5/SD)*sizes
            ztd = (-6*sin(xtd*TD*pi/50)/TD-1.5*sin((xtd*TD*pi/25)+PROFT)/TD+7.5/TD)*sizet
            z = zfd + max(zsd, ztd) + DPOSIT
            z = max(z, ((7.5/FD)+(7.5/SD)+(7.5/TD))*self.ELVMIN+DPOSIT)
        elif self.TYPE == 6:
            zfd = (-6*sin(xfd*FD*pi/50)/FD-1.5*sin((xfd*FD*pi/25)+PROFF)/FD)
            shape = ((7.5/FD)-zfd)/(15/FD)
            zfd = zfd*sizef
            zsd = (-6*sin(xsd*SD*pi/50)/SD-1.5*sin((xsd*SD*pi/25)+PROFS)/SD)*shape*sizes
            ztd = (-6*sin(xtd*TD*pi/50)/TD-1.5*sin((xtd*TD*pi/25)+PROFT)/TD)*shape*sizet
            z = zfd + zsd + ztd + DPOSIT
            z = max(z, ((7.5 / FD) + (7.5 / SD) + (7.5 / TD)) * self.ELVMIN + DPOSIT)
        else:
            pass

        # Determine high and low points on bedform surface, and set initial values
        # of elevation arrays (zcont and ZBED).
        if run == 0:
            zmin = min(min(z))
            zmax = max(max(z))
            zref = zmin + (zmax-zmin) * self.ZHORIZ
            if zref < -30.0:
                zref = -30.0
            ZBED[GridSize/EdgeGridSpace:-EdgeGridSpace:1] = z[:,1]
            ZBED[GridSize/EdgeGridSpace:2*GridSize/EdgeGridSpace-1] = z[1,:]
            ZBED[2*GridSize/EdgeGridSpace:3*GridSize/EdgeGridSpace-1] = z[:,100]
            ZBED[4*GridSize/EdgeGridSpace-1:-EdgeGridSpace:3*GridSize/EdgeGridSpace] = z[100,:]
            ZBED = min(ZBED,zref)
            ZBED = max(ZBED,-30)
        # Define ZBED
        if np.size(z, 1) > 1:
            ZBED[GridSize/EdgeGridSpace:-EdgeGridSpace:1] = \
                min(ZBED[GridSize/EdgeGridSpace:-EdgeGridSpace:1],z[:,1].T)
            ZBED[GridSize/EdgeGridSpace:2*GridSize/EdgeGridSpace-1] = \
                min(ZBED[GridSize/EdgeGridSpace:2*GridSize/EdgeGridSpace-1],z[1,:])
            ZBED[2*GridSize/EdgeGridSpace:3*GridSize/EdgeGridSpace-1] = \
                min(ZBED[2*GridSize/EdgeGridSpace:3*GridSize/EdgeGridSpace-1],z[:,100].T)
        else:
            ZBED = min(z, ZBED)
            ZBED = max(ZBED,-30)
        return ZBED