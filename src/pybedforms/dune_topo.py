import numpy as np
from numpy import pi, sin, cos, arctan2, tan
import numba as nb
from tqdm import tqdm
from mayavi import mlab
import matplotlib.pyplot as plt
from skimage import measure
import scipy

def make_grid_and_edges(center = 50, size = 100, grid_spacing = 1, edge_grid_spacing = 1):
    x, y = np.meshgrid(
        np.arange(1 - center, size - center + grid_spacing, grid_spacing),
        np.arange(1 - center, size - center + grid_spacing, grid_spacing)
    )
    num_cells = size // edge_grid_spacing
    #todo: still off by 1 relative to matlab but this is much cleaner
    XEdge = np.concatenate((
        np.ones(num_cells),
        np.arange(1, size + edge_grid_spacing, edge_grid_spacing),
        np.repeat(size, num_cells),
        np.arange(1, size + edge_grid_spacing, edge_grid_spacing)[::-1]
    )) - center
    YEdge = np.concatenate((
        np.arange(1, size + edge_grid_spacing, edge_grid_spacing)[::-1],
        np.ones(num_cells),
        np.arange(1, size + edge_grid_spacing, edge_grid_spacing),
        np.repeat(size, num_cells),
    )) - center
    return x, y, XEdge, YEdge

class DuneTopo(object):
    def __init__(self, SPCNGF=50.0, PHASEF=0.0, SMTRYF=1.0, SMCHGF=0.0, SMPRDF=1.0,
                 SMFAZF=0.0,HTRTOF=1.0,HTCHGF=0.0,HTPRDF=1.0,HTFAZF=0.0,SNSPF1=0.0,
                 SNMGF1=0.0,SNFZF1=0.0,SNVLF1=0.0,SNSPF2=0.0,SNMGF2=0.0,SNFZF2=0.0,
                 SNVLF2=0.0,TRENDF=90.0,VELOCF=1.0,VLCHGF=0.0,VLPRDF=1.0,VLFAZF=0.0,
                 SPCNGS=0.0,PHASES=0.0,SMTRYS=0.0,SMCHGS=0.0,SMPRDS=1.0,SMFAZS=0.0,
                 HTRTOS=0.0,HTCHGS=0.0,HTPRDS=1.0,HTFAZS=0.0,SNSPS1=0.0,SNMGS1=0.0,
                 SNFZS1=0.0,SNVLS1=0.0,SNSPS2=0.0,SNMGS2=0.0,SNFZS2=0.0,SNVLS2=0.0,
                 TRENDS=0.0,VELOCS=0.0,VLCHGS=1.0,VLPRDS=1.0,VLFAZS=0.0,SPCNGT=0.0,
                 PHASET=0.0,SMTRYT=0.0,SMCHGT=0.0,SMPRDT=1.0,SMFAZT=0.0,HTRTOT=0.0,
                 HTCHGT=0.0,HTPRDT=1.0,HTFAZT=0.0,SNSPT1=0.0,SNMGT1=0.0,SNFZT1=0.0,
                 SNVLT1=0.0,SNSPT2=0.0,SNMGT2=0.0,SNFZT2=0.0,SNVLT2=0.0,TRENDT=0.0,
                 VELOCT=0.0,VLCHGT=0.0,VLPRDT=1.0,VLFAZT=0.0,TYPE=1,CHOICE=0,
                 ELVMIN=-1.0,DEPRAT=.08,DEPCHG=0.0,DEPPRD=1.0,DEPFAZ=0.0,NBEDSH=500,
                 INTXBD=4,FRMNUM=1,CAPLOG='true',FILOG='true',FRMLOG='true',IGRDSP=1,
                 ZHORIZ=1.0):
        # default parameter values are from 'Figure 5' of Rubin and Carter
        self.TOPO = None
        self.SPCNGF = SPCNGF # (3)  WAVELENGTH OF BEDFORMS IN FIRST SET
        self.PHASEF = PHASEF # (4)  Bedform phase (controls placement within block diagram)
        self.SMTRYF = SMTRYF # (5)  Mean asymmetry
        self.SMCHGF = SMCHGF # (6)  Amplitude of asymmetry cycle
        self.SMPRDF = SMPRDF # (7)  Period of asymmetry cycle
        self.SMFAZF = SMFAZF # (8)  Phase of asymmetry cycle
        self.HTRTOF = HTRTOF # (9)  Mean steepness
        self.HTCHGF = HTCHGF # (10) Amplitude of steepness cycle
        self.HTPRDF = HTPRDF # (11) Period of steepness cycle
        self.HTFAZF = HTFAZF # (12) Phase of steepness cycle
        self.SNSPF1 = SNSPF1 # (13) Wavelength (alon g -crest) of first set of pla n -form sinuosities
        self.SNMGF1 = SNMGF1 # (14) Amplitude (measured in pla n -form) of first set of sinuosities
        self.SNFZF1 = SNFZF1 # (15) Phase of first set of sinuosities
        self.SNVLF1 = SNVLF1 # (16) Migration speed (alon g -crest) of first set of sinuosities
        self.SNSPF2 = SNSPF2 # (17) Wavelength (alon g -crest) of second set of pla n -form sinuosities
        self.SNMGF2 = SNMGF2 # (18) Amplitude (measured in pla n -form) of second set of sinuosities
        self.SNFZF2 = SNFZF2 # (19) Phase of second set of sinuosities
        self.SNVLF2 = SNVLF2 # (20) Migration speed (alon g -crest) of second set of sinuosities
        self.TRENDF = TRENDF # (21) Migration direction of bedform
        self.VELOCF = VELOCF # (22) Mean migration speed of bedform
        self.VLCHGF = VLCHGF # (23) Amplitude of speed cycle
        self.VLPRDF = VLPRDF # (24) Period of speed cycle
        self.VLFAZF = VLFAZF # (25) Phase of speed cycle
        self.SPCNGS = SPCNGS # (26) WAVELENGTH OF BEDFORMS IN SECOND SET
        self.PHASES = PHASES # (27) Bedform phase (controls placement within block diagram)
        self.SMTRYS = SMTRYS # (28) Mean asymmetry
        self.SMCHGS = SMCHGS # (29) Amplitude of asymmetry cycle
        self.SMPRDS = SMPRDS # (30) Period of asymmetry cycle
        self.SMFAZS = SMFAZS # (31) Phase of asymmetry cycle
        self.HTRTOS = HTRTOS # (32) Mean steepness
        self.HTCHGS = HTCHGS # (33) Amplitude of steepness cycle
        self.HTPRDS = HTPRDS # (34) Period of steepness cycle
        self.HTFAZS = HTFAZS # (35) Phase of steepness cycle
        self.SNSPS1 = SNSPS1 # (36) Wavelength (alon g -crest) of first set of pla n -form sinuosities
        self.SNMGS1 = SNMGS1 # (37) Amplitude (measured in pla n -form) of first set of sinuosities
        self.SNFZS1 = SNFZS1 # (38) Phase of first set of sinuosities
        self.SNVLS1 = SNVLS1 # (39) Migration speed (alon g -crest) of first set of sinuosities
        self.SNSPS2 = SNSPS2 # (40) Wavelength (alon g -crest) of second set of pla n -form sinuosities
        self.SNMGS2 = SNMGS2 # (41) Amplitude (measured in pla n -form) of second set of sinuosities
        self.SNFZS2 = SNFZS2 # (42) Phase of second set of sinuosities
        self.SNVLS2 = SNVLS2 # (43) Migration speed (alon g -crest) of second set of sinuosities
        self.TRENDS = TRENDS # (44) Migration direction of bedform
        self.VELOCS = VELOCS # (45) Mean migration speed of bedform
        self.VLCHGS = VLCHGS # (46) Amplitude of speed cycle
        self.VLPRDS = VLPRDS # (47) Period of speed cycle
        self.VLFAZS = VLFAZS # (48) Phase of speed cycle
        self.SPCNGT = SPCNGT # (49) WAVELENGTH OF BEDFORMS IN THIRD SET
        self.PHASET = PHASET # (50) Bedform phase (controls placement within block diagram)
        self.SMTRYT = SMTRYT # (51) Mean asymmetry
        self.SMCHGT = SMCHGT # (52) Amplitude of asymmetry cycle
        self.SMPRDT = SMPRDT # (53) Period of asymmetry cycle
        self.SMFAZT = SMFAZT # (54) Phase of asymmetry cycle
        self.HTRTOT = HTRTOT # (55) Mean steepness
        self.HTCHGT = HTCHGT # (56) Amplitude of steepness cycle
        self.HTPRDT = HTPRDT # (57) Period of steepness cycle
        self.HTFAZT = HTFAZT # (58) Phase of steepness cycle
        self.SNSPT1 = SNSPT1 # (59) Wavelength (alon g -crest) of first set of pla n -form sinuosities
        self.SNMGT1 = SNMGT1 # (60) Amplitude (measured in pla n -form) of first set of sinuosities
        self.SNFZT1 = SNFZT1 # (61) Phase of first set of sinuosities
        self.SNVLT1 = SNVLT1 # (62) Migration speed (alon g -crest) of first set of sinuosities
        self.SNSPT2 = SNSPT2 # (63) Wavelength (alon g -crest) of second set of pla n -form sinuosities
        self.SNMGT2 = SNMGT2 # (64) Amplitude (measured in pla n -form) of second set of sinuosities
        self.SNFZT2 = SNFZT2 # (65) Phase of second set of sinuosities
        self.SNVLT2 = SNVLT2 # (66) Migration speed (alon g -crest) of second set of sinuosities
        self.TRENDT = TRENDT # (67) Migration direction of bedform
        self.VELOCT = VELOCT # (68) Mean migration speed of bedform
        self.VLCHGT = VLCHGT # (69) Amplitude of speed cycle
        self.VLPRDT = VLPRDT # (70) Period of speed cycle
        self.VLFAZT = VLFAZT # (71) Phase of speed cycle
        self.TYPE   = TYPE   # (72) Type of superpositioning (INTEGER)
        self.CHOICE = CHOICE # (73) Rotation option (INTEGER
        self.ELVMIN = ELVMIN # (74) Elevation of interdune flats
        self.DEPRAT = DEPRAT # (75) Rate of deposition
        self.DEPCHG = DEPCHG # (76) Amplitude of cycle in rate of deposition
        self.DEPPRD = DEPPRD # (77) Period of cycle in rate of deposition
        self.DEPFAZ = DEPFAZ # (78) Phase of cycle in rate of deposition
        self.NBEDSH = NBEDSH # (79) Time from t=0 to beginning of depositional episode (INTEGER)
        self.INTXBD = INTXBD # (80) Interval between drawing of crossbeds (INTEGER)
        self.FRMNUM = FRMNUM # (81) Time from t=0 to end of depositional episode (INTEGER)
        self.CAPLOG = CAPLOG # (82) Print caption? (LOGICAL VARIABLE)
        self.FILOG  = FILOG  # (83) Print name of input file? (LOGICAL VARIABLE)
        self.FRMLOG = FRMLOG # (84) Print time at end of depositional episode? (LOGICAL VARIABLE)
        self.IGRDSP = IGRDSP # (85) Precision? (low number) or speed? (high number) (1 2 4 5 or 10)
        self.ZHORIZ = ZHORIZ # (86) Elevation of horizontal section
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

        # Make the grids from dune init
        self.CenterShift = 50
        self.GridSize = 100
        self.SurfGridSpace = 1
        self.EdgeGridSpace = 1
        self.x, self.y, self.XEdge, self.YEdge = make_grid_and_edges(self.CenterShift, self.GridSize, self.SurfGridSpace, self.EdgeGridSpace)
        """
        % Calculate the migration direction of scour pits formed by intersecting
        % bedform troughs of the first two sets of bedforms.  When specified in 
        % the input paramaters, this calculation is used to rotate the bedforms 
        % such that the sides of the block diagram are normal and parallel to the 
        % axes of trough-shaped sets of cross-bedding..
        """
        if (tan(self.TRENDF * 2 * pi / 360) != tan(self.TRENDS * 2 * pi / 360)) & (self.SPCNGF > 0) & (self.SPCNGS > 0):
            self.ANGLEA = arctan2(
                self.VELOCF, tan((90 + self.TRENDF - self.TRENDS) * 2 * pi / 360) * self.VELOCF - self.VELOCS / sin(
                    (self.TRENDS - self.TRENDF) * 2 * pi / 360)
            ) * 360 / (2 * pi) - 90
        if self.CHOICE != 0:
            self.ANGLEB = self.TRENDS - self.TRENDF
            self.ANGLEC = self.TRENDT - self.TRENDF
            self.TRENDF -= self.ANGLEA + (self.CHOICE - 1) * 90
            self.TRENDS = self.TRENDF + self.ANGLEB
            if self.SPCNGT != 0:
                self.TRENDT = self.TRENDF + self.ANGLEC

    def run(self, TIME): #, run): #, GridSize, EdgeGridSpace):
        """
        main run method, below will be split into other functions as needed
        :return:
        """
        # Calculate deposition.
        DPOSIT = (self.DEPRAT * TIME) + self.DEPCHG * (-self.DEPPRD / (2 * pi)) * cos(
            ((self.DEPFAZ - 90) * 2 * pi / 360) + (TIME * 2 * pi / self.DEPPRD))

        # Calculate bedform asymmetry.
        PROFF = ((1 - self.SMTRYF) + self.SMCHGF * sin(((-TIME / self.SMPRDF) - self.SMFAZF / 360) * pi * 2.0)) * pi / 2.0
        PROFS = ((1 - self.SMTRYS) + self.SMCHGS * sin(((-TIME / self.SMPRDS) - self.SMFAZS / 360) * pi * 2.0)) * pi / 2.0
        PROFT = ((1 - self.SMTRYT) + self.SMCHGT * sin(((-TIME / self.SMPRDT) - self.SMFAZT / 360) * pi * 2.0)) * pi / 2.0

        # Calculate bedform sizes.
        sizef = self.HTRTOF + (self.HTCHGF * sin(((TIME / self.HTPRDF) * 2 * pi) + ((self.HTFAZF - 90) * 2 * pi / 360.0)))
        sizes = self.HTRTOS + (self.HTCHGS * sin(((TIME / self.HTPRDS) * 2 * pi) + ((self.HTFAZS - 90) * 2 * pi / 360.0)))
        sizet = self.HTRTOT + (self.HTCHGT * sin(((TIME / self.HTPRDT) * 2 * pi) + ((self.HTFAZT - 90) * 2 * pi / 360.0)))

        # Calculate bedform locations.
        dispfd = self.VELOCF * TIME + self.VLCHGF * (-self.VLPRDF / (2 * pi)) * cos(
            ((self.VLFAZF - 90) * 2 * pi / 360) + (TIME * 2 * pi / self.VLPRDF))
        dispsd = self.VELOCS * TIME + self.VLCHGS * (-self.VLPRDS / (2 * pi)) * cos(
            ((self.VLFAZS - 90) * 2 * pi / 360) + (TIME * 2 * pi / self.VLPRDS))
        disptd = self.VELOCT * TIME + self.VLCHGT * (-self.VLPRDT / (2 * pi)) * cos(
            ((self.VLFAZT - 90) * 2 * pi / 360) + (TIME * 2 * pi / self.VLPRDT))

        yfd = self.x * sin(self.TRENDF * 2 * pi / 360) + self.y * cos(self.TRENDF * 2 * pi / 360)
        ysd = self.x * sin(self.TRENDS * 2 * pi / 360) + self.y * cos(self.TRENDS * 2 * pi / 360)
        ytd = self.x * sin(self.TRENDT * 2 * pi / 360) + self.y * cos(self.TRENDT * 2 * pi / 360)
        xfd = (self.x * cos(self.TRENDF * 2 * pi / 360) - self.y * sin(self.TRENDF * 2 * pi / 360) - dispfd - self.SPCNGF * self.PHASEF / 360) - \
              (self.SNMGF1 * sin(
                  (yfd * 2 * pi / self.SNSPF1) + ((self.SNFZF1 * 2 * pi / 360) + (TIME * self.SNVLF1 * 2 * pi / self.SNSPF1)))) - \
              (self.SNMGF2 * sin(
                  (yfd * 2 * pi / self.SNSPF2) + ((self.SNFZF2 * 2 * pi / 360) + (TIME * self.SNVLF2 * 2 * pi / self.SNSPF2))))
        xsd = (self.x * cos(self.TRENDS * 2 * pi / 360) - self.y * sin(self.TRENDS * 2 * pi / 360) - dispsd - self.SPCNGS * self.PHASES / 360) - \
              (self.SNMGS1 * sin(
                  (ysd * 2 * pi / self.SNSPS1) + ((self.SNFZS1 * 2 * pi / 360) + (TIME * self.SNVLS1 * 2 * pi / self.SNSPS1)))) - \
              (self.SNMGS2 * sin(
                  (ysd * 2 * pi / self.SNSPS2) + ((self.SNFZS2 * 2 * pi / 360) + (TIME * self.SNVLS2 * 2 * pi / self.SNSPS2))))
        xtd = (self.x * cos(self.TRENDT * 2 * pi / 360) - self.y * sin(self.TRENDT * 2 * pi / 360) - disptd - self.SPCNGT * self.PHASET / 360) - \
              (self.SNMGT1 * sin(
                  (ytd * 2 * pi / self.SNSPT1) + ((self.SNFZT1 * 2 * pi / 360) + (TIME * self.SNVLT1 * 2 * pi / self.SNSPT1)))) - \
              (self.SNMGT2 * sin(
                  (ytd * 2 * pi / self.SNSPT2) + ((self.SNFZT2 * 2 * pi / 360) + (TIME * self.SNVLT2 * 2 * pi / self.SNSPT2))))
        # Superimpose sets of dunes
        FD = self.FD
        SD = self.SD
        TD = self.TD
        z = None
        if self.TYPE == 1:
            zfd = (-6 * sin(xfd * FD * pi / 50) / FD - 1.5 * sin((xfd * FD * pi / 25) + PROFF) / FD)
            shape = 1
            zfd = zfd * sizef
            zsd = (-6 * sin(xsd * SD * pi / 50) / SD - 1.5 * sin((xsd * SD * pi / 25) + PROFS) / SD) * shape * sizes
            ztd = (-6 * sin(xtd * TD * pi / 50) / TD - 1.5 * sin((xtd * TD * pi / 25) + PROFT) / TD) * shape * sizet
            z = zfd + zsd + ztd + DPOSIT
            z = np.maximum(z, ((7.5 / FD) + (7.5 / SD) + (7.5 / TD)) * self.ELVMIN + DPOSIT)
        elif self.TYPE == 2:
            zfd = (-6 * sin(xfd * FD * pi / 50) / FD - 1.5 * sin((xfd * FD * pi / 25) + PROFF) / FD)
            shape = ((7.5 / FD) - zfd) / (15 / FD)
            zfd = zfd * sizef
            zsd = (-6 * sin(xsd * SD * pi / 50) / SD - 1.5 * sin((xsd * SD * pi / 25) + PROFS) / SD) * shape * sizes
            ztd = (-6 * sin(xtd * TD * pi / 50) / TD - 1.5 * sin((xtd * TD * pi / 25) + PROFT) / TD) * shape * sizet
            z = zfd + zsd + ztd + DPOSIT
            z = np.maximum(z, ((7.5 / FD) + (7.5 / SD) + (7.5 / TD)) * self.ELVMIN + DPOSIT)
        elif self.TYPE == 3:
            BD = 100 / max(self.SPCNGF, self.SPCNGS, self.SPCNGT)
            zfd = (-6 * sin(xfd * FD * pi / 50) / FD - 1.5 * sin((xfd * FD * pi / 25) + PROFF) / FD + 7.5 / FD) * sizef
            zsd = (-6 * sin(xsd * SD * pi / 50) / SD - 1.5 * sin((xsd * SD * pi / 25) + PROFS) / SD + 7.5 / SD) * sizes
            ztd = (-6 * sin(xtd * TD * pi / 50) / TD - 1.5 * sin((xtd * TD * pi / 25) + PROFT) / TD + 7.5 / TD) * sizet
            z = np.maximum(zfd, zsd)
            z = np.maximum(z, ztd)
            z = z + DPOSIT
            z = np.maximum(z, (7.5 * (1 + self.ELVMIN) / BD + DPOSIT))
        elif self.TYPE == 4:
            zfd = (-6 * sin(xfd * FD * pi / 50) / FD - 1.5 * sin((xfd * FD * pi / 25) + PROFF) / FD)
            shape = 1 - ((7.5 / FD) - zfd) / (15 / FD)
            zfd = zfd * sizef
            zsd = (-6 * sin(xsd * SD * pi / 50) / SD - 1.5 * sin((xsd * SD * pi / 25) + PROFS) / SD) * shape * sizes
            ztd = (-6 * sin(xtd * TD * pi / 50) / TD - 1.5 * sin((xtd * TD * pi / 25) + PROFT) / TD) * shape * sizet
            z = zfd + zsd + ztd + DPOSIT
            z = np.maximum(z, ((7.5 / FD) + (7.5 / SD) + (7.5 / TD)) * self.ELVMIN + DPOSIT)
        elif self.TYPE == 5:
            zfd = (-6 * sin(xfd * FD * pi / 50) / FD - 1.5 * sin((xfd * FD * pi / 25) + PROFF) / FD) * sizef
            zsd = (-6 * sin(xsd * SD * pi / 50) / SD - 1.5 * sin((xsd * SD * pi / 25) + PROFS) / SD + 7.5 / SD) * sizes
            ztd = (-6 * sin(xtd * TD * pi / 50) / TD - 1.5 * sin((xtd * TD * pi / 25) + PROFT) / TD + 7.5 / TD) * sizet
            z = zfd + np.maximum(zsd, ztd) + DPOSIT
            z = np.maximum(z, ((7.5 / FD) + (7.5 / SD) + (7.5 / TD)) * self.ELVMIN + DPOSIT)
        elif self.TYPE == 6:
            zfd = (-6 * sin(xfd * FD * pi / 50) / FD - 1.5 * sin((xfd * FD * pi / 25) + PROFF) / FD)
            shape = ((7.5 / FD) - zfd) / (15 / FD)
            zfd = zfd * sizef
            zsd = (-6 * sin(xsd * SD * pi / 50) / SD - 1.5 * sin((xsd * SD * pi / 25) + PROFS) / SD) * shape * sizes
            ztd = (-6 * sin(xtd * TD * pi / 50) / TD - 1.5 * sin((xtd * TD * pi / 25) + PROFT) / TD) * shape * sizet
            z = zfd + zsd + ztd + DPOSIT
            z = np.maximum(z, ((7.5 / FD) + (7.5 / SD) + (7.5 / TD)) * self.ELVMIN + DPOSIT)
        else:
            pass

        if TIME == 0:
            self.TOPO = z
        else:
            self.TOPO = np.dstack((self.TOPO, z))

        # I am storing all z surfaces as 'TOPO' (so that we have a complete stratigraphic model), so no need for 'ZBED'
        # Determine high and low points on bedform surface, and set initial values
        # of elevation arrays (zcont and ZBED).
        # n_gridcells = int(GridSize/EdgeGridSpace)  
        # if run == 0:
        #     zmin = np.min(z)
        #     zmax = np.max(z)
        #     zref = zmin + (zmax - zmin) * self.ZHORIZ
        #     if zref < -30.0:
        #         zref = -30.0
        #     ZBED[n_gridcells : 0 : -EdgeGridSpace] = z[:, 0]
        #     ZBED[n_gridcells : 2 * n_gridcells] = z[0, :]
        #     ZBED[2 * n_gridcells : 3 * n_gridcells] = z[:, -1]
        #     ZBED[4 * n_gridcells : 3 * n_gridcells-1 : -EdgeGridSpace] = z[-1, :]
        #     ZBED = np.minimum(ZBED, zref)
        #     ZBED = np.maximum(ZBED, -30)
        # # Define ZBED
        # if np.size(z, 0) > 1:
        #     ZBED[n_gridcells : 0 : -EdgeGridSpace] = \
        #         np.minimum(ZBED[n_gridcells : 0 : -EdgeGridSpace], z[:, 0].T)
        #     ZBED[n_gridcells : 2 * n_gridcells] = \
        #         np.minimum(ZBED[n_gridcells : 2 * n_gridcells], z[0, :])
        #     ZBED[2 * n_gridcells : 3 * n_gridcells] = \
        #         np.minimum(ZBED[2 * n_gridcells : 3 * n_gridcells], z[:, -1].T)
        # else:
        #     ZBED = np.minimum(z, ZBED)
        #     ZBED = np.maximum(ZBED, -30)
        # return ZBED

def plot_3d(topo, bottom_z, top_z, time_step, color=(0.88627451, 0.79215686, 0.4627451), scale=1, ve=1, dx=1, 
                                                    line_thickness=0.05, contour_switch=False, new_figure=False):
    """function for plotting a set of bedforms in 3D"""

    if new_figure:
        mlab.figure(bgcolor=(1,1,1))
    else:
        mlab.clf()
        
    strat = np.minimum.accumulate(topo[:, :, ::-1], axis=2)[:, :, ::-1] # convert topography to stratigraphy
    strat2 = strat.copy()
    strat2[strat<bottom_z] = bottom_z
    strat2[strat>top_z] = top_z
    
    r,c,ts = np.shape(strat2)
    z = scale*strat2[:,:,ts-1]
    z1 = strat2[:,:,-1]

    X1 = scale*(np.linspace(0,r-1,r)*dx)
    Y1 = scale*(np.linspace(0,c-1,c)*dx)
    mlab.surf(X1,Y1,z,warp_scale=ve,color=color)

    if contour_switch:
        contours = list(np.arange(vmin,vmax,ci*scale)) # list of contour values
        mlab.contour_surf(X1,Y1,z,contours=contours,warp_scale=ve,color=(0,0,0),line_width=1.0)

    # updip side:
    vertices, triangles = create_section(z1[:,0],dx,bottom_z) 
    x = scale*(vertices[:,0])
    y = scale*(np.zeros(np.shape(vertices[:,0])))
    z = scale*ve*vertices[:,1]
    mlab.triangular_mesh(x,y,z,triangles,color=color)

    # downdip side:
    vertices, triangles = create_section(z1[:,-1],dx,bottom_z) 
    x = scale*(vertices[:,0])
    y = scale*((c-1)*dx*np.ones(np.shape(vertices[:,0])))
    z = scale*ve*vertices[:,1]
    mlab.triangular_mesh(x,y,z,triangles,color=color)

    # left edge (looking downdip):
    vertices, triangles = create_section(z1[0,:],dx,bottom_z) 
    x = scale*(np.zeros(np.shape(vertices[:,0])))
    y = scale*(vertices[:,0])
    z = scale*ve*vertices[:,1]
    mlab.triangular_mesh(x,y,z,triangles,color=color)

    # right edge (looking downdip):
    vertices, triangles = create_section(z1[-1,:],dx,bottom_z) 
    x = scale*((r-1)*dx*np.ones(np.shape(vertices[:,0])))
    y = scale*(vertices[:,0])
    z = scale*ve*vertices[:,1]
    mlab.triangular_mesh(x,y,z,triangles,color=color)

    # bottom face of block:
    vertices = dx*np.array([[0,0],[r-1,0],[r-1,c-1],[0,c-1]])
    triangles = [[0,1,3],[1,3,2]]
    x = scale*(vertices[:,0])
    y = scale*(vertices[:,1])
    z = scale*bottom_z*np.ones(np.shape(vertices[:,0]))
    mlab.triangular_mesh(x,y,ve*z,triangles,color=color)

    t_steps = np.hstack((np.arange(0, ts-1, time_step), ts-2))
    for layer_n in tqdm(t_steps): # main loop
        top = strat2[:,0,layer_n+1]  # updip side
        base = strat2[:,0,layer_n]
        X1 = scale*(dx*np.arange(0,r))
        Y1 = scale*(np.zeros(np.shape(base)))
        Z1 = ve*scale*base
        mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)
        if layer_n == ts-2:
            Z1 = ve*scale*top
            mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)

        top = strat2[:,-1,layer_n+1]  # downdip side
        base = strat2[:,-1,layer_n]
        X1 = scale*(dx*np.arange(0,r))
        Y1 = scale*(dx*(c-1)*np.ones(np.shape(base)))
        Z1 = ve*scale*base
        mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)
        if layer_n == ts-2:
            Z1 = ve*scale*top
            mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)

        top = strat2[0,:,layer_n+1]  # left edge (looking downdip)
        base = strat2[0,:,layer_n]
        X1 = scale*(np.zeros(np.shape(base)))
        Y1 = scale*(dx*np.arange(0,c))
        Z1 = ve*scale*base
        mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)
        if layer_n == ts-2:
            Z1 = ve*scale*top
            mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)

        top = strat2[-1,:,layer_n+1] # right edge (looking downdip)
        base = strat2[-1,:,layer_n]
        X1 = scale*(dx*(r-1)*np.ones(np.shape(base)))
        Y1 = scale*(dx*np.arange(0,c))
        Z1 = ve*scale*base
        mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)
        if layer_n == ts-2:
            Z1 = ve*scale*top
            mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)

        if (bottom_z > np.min(strat[:,:,layer_n])) and (bottom_z < np.max(strat[:,:,layer_n])):
            contours = measure.find_contours(strat[:,:,layer_n], bottom_z)
            for i in range(0,len(contours)):
                x1 = contours[i][:,0]
                y1 = contours[i][:,1]
                z1 = np.ones(np.shape(x1)) * bottom_z
                mlab.plot3d(x1,y1,z1,color=(0,0,0),tube_radius=line_thickness)

        if (top_z > np.min(strat[:,:,layer_n])) and (top_z < np.max(strat[:,:,layer_n])):
            contours = measure.find_contours(strat[:,:,layer_n], top_z)
            for i in range(0,len(contours)):
                x1 = contours[i][:,0]
                y1 = contours[i][:,1]
                z1 = np.ones(np.shape(x1)) * top_z
                mlab.plot3d(x1,y1,z1,color=(0,0,0),tube_radius=line_thickness)

def create_section(profile, dx, bottom):
    """function for creating a cross section from a top surface
    inputs:
    profile - elevation data for top surface
    dx - gridcell size
    bottom - elevation value for the bottom of the block
    returns:
    vertices - coordinates of vertices
    triangles - indices of the 'vertices' array that from triangles (for triangular mesh)
    """
    x1 = dx*np.linspace(0, len(profile)-1, len(profile))
    x = np.hstack((x1, x1[::-1]))
    y = np.hstack((profile, bottom*np.ones(np.shape(x1))))
    vertices = np.vstack((x, y)).T
    n = len(x)
    triangles = []
    for i in range(0,int((n-1)/2)):
        triangles.append([i,i+1,n-i-1])
        triangles.append([i+1,n-i-1,n-i-2])
    return vertices, triangles

def extract_core(topo, scale = 1, ve = 1, x0 = 50, y0 = 50, dx = 1, radius = 3, num = 50, bottom = 0, time_step = 4):
    """function for extracting a "core" from the model"""
    strat = np.minimum.accumulate(topo[:, :, ::-1], axis=2)[:, :, ::-1] # convert topography to stratigraphy
    r, c, ts = np.shape(strat)
    strat2 = strat.copy()
    strat2[strat<bottom] = bottom
    X1 = x0 + np.cos(2*pi/num*np.arange(num))*radius
    X1 = np.hstack((X1, X1[0]))
    Y1 = y0 + np.sin(2*pi/num*np.arange(num))*radius
    Y1 = np.hstack((Y1, Y1[0]))
    top = scipy.ndimage.map_coordinates(strat2[:,:,-1], np.vstack((Y1, X1)))
    vertices, triangles = create_section(top, dx, bottom) 
    color = (0.886, 0.792, 0.463) # color for plotting basal part of panel
    mlab.triangular_mesh(scale*np.hstack((dx*X1,dx*X1[::-1])),
                        scale*np.hstack((dx*Y1,dx*Y1[::-1])),scale*ve*vertices[:,1],triangles,color=color)
    t_steps = np.hstack((np.arange(0, ts, time_step), ts-1))
    for layer_n in tqdm(t_steps): # main loop
        Z1 = scipy.ndimage.map_coordinates(strat2[:,:,layer_n], np.vstack((Y1, X1)))
        mlab.plot3d(X1, Y1, Z1, color=(0,0,0), tube_radius=0.05)