#! /usr/bin/env python
#
# Time-stamp: <2016-02-12 16:03:47 hamana>
#
#------------------------
# sip2tpv.py based on SIPtoTPV.py
#------------------------
#
# usage: sip2tpv.py input[wcs.fits] output[filename.head]
#
# Copyright: Junya Sakurai, Takashi Hamana, Satoshi Miyazaki (National Astronomical Observatory Japan)
# 
# Note: This is a free software. 
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Importrant note:
#       sip2tpv.py can transform 9th order SIP to 9th order TPV,
#       but 9th order TPV is NOT officially supported.
#       Most of standard softwares cannot treat 9th order TPV.
#
#       If you publish a paper using this software, 
#       cite Hamana et al. (2015) (Publ Astron Soc Jpn (June 2015) 67 (3), 34, doi:10.1093/pasj/psv034)
# 

import sys
import logging
import argparse
from astropy.io import fits
import numpy as np
import os.path
from collections import OrderedDict

def sip2tpv(inhdr):
    """Convert the SIP coefficients contained in dictionary inhdr to
    TPV coefficients. Return the TPV coefficients as an OrderedDict.
    """
    # check SIP order 
    Aorder=int(inhdr['A_ORDER'])
    if Aorder>7:
        logging.warn(f'A_ORDER = {Aorder}')
        logging.warn('Be sure, this is out of the specification of TPV')

    if Aorder>9:
        logging.error(f'A_ORDER = {Aorder}')
        logging.error('Bad A_ORDER, abort')
        raise ValueError(f'Bad A_ORDER = {Aorder}')


    Border=int(inhdr['B_ORDER'])
    if Border>7:
        logging.warn(f'B_ORDER = {Border}')
        logging.warn('Be sure, this is out of the specification of TPV')

    if Border>9:
        logging.error('B_ORDER = {Border}')
        logging.error('Bad B_ORDER, abort')
        raise ValueError(f'Bad B_ORDER = {Border}')

    outhdr = OrderedDict()

    #WCS parameters#
    if 'EQUINOX' in inhdr:
        outhdr['EQUINOX'] = inhdr['EQUINOX']
    outhdr['RADESYS'] = inhdr['RADESYS']
    outhdr['CRPIX1']  = inhdr['CRPIX1']
    outhdr['CRPIX2']  = inhdr['CRPIX2']
    outhdr['CRVAL1']  = inhdr['CRVAL1']
    outhdr['CRVAL2']  = inhdr['CRVAL2']
    if 'CUNIT1' in inhdr:
        outhdr['CUNIT1']  = inhdr['CUNIT1']
    if 'CUNIT2' in inhdr:
        outhdr['CUNIT2']  = inhdr['CUNIT2']
    outhdr['CTYPE1']  = 'RA---TPV'
    outhdr['CTYPE2']  = 'DEC--TPV'

    #CD matrix#
    C11=inhdr['CD1_1']
    C12=inhdr['CD1_2']
    C21=inhdr['CD2_1']
    C22=inhdr['CD2_2']

    outhdr['CD1_1'] = C11
    outhdr['CD1_2'] = C12
    outhdr['CD2_1'] = C21
    outhdr['CD2_2'] = C22

    #coefficient A#
    if Aorder>=2:
        A20=inhdr['A_2_0']
        A11=inhdr['A_1_1']
        A02=inhdr['A_0_2']

    if Aorder>=3:
        A30=inhdr['A_3_0']
        A21=inhdr['A_2_1']
        A12=inhdr['A_1_2']
        A03=inhdr['A_0_3']

    if Aorder>=4:
        A40=inhdr['A_4_0']
        A31=inhdr['A_3_1']
        A22=inhdr['A_2_2']
        A13=inhdr['A_1_3']
        A04=inhdr['A_0_4']

    if Aorder>=5:
        A50=inhdr['A_5_0']
        A41=inhdr['A_4_1']
        A32=inhdr['A_3_2']
        A23=inhdr['A_2_3']
        A14=inhdr['A_1_4']
        A05=inhdr['A_0_5']

    if Aorder>=6:
        A60=inhdr['A_6_0']
        A51=inhdr['A_5_1']
        A42=inhdr['A_4_2']
        A33=inhdr['A_3_3']
        A24=inhdr['A_2_4']
        A15=inhdr['A_1_5']
        A06=inhdr['A_0_6']

    if Aorder>=7:
        A70=inhdr['A_7_0']
        A61=inhdr['A_6_1']
        A52=inhdr['A_5_2']
        A43=inhdr['A_4_3']
        A34=inhdr['A_3_4']
        A25=inhdr['A_2_5']
        A16=inhdr['A_1_6']
        A07=inhdr['A_0_7']

    if Aorder>=8:
        A80=inhdr['A_8_0']
        A71=inhdr['A_7_1']
        A62=inhdr['A_6_2']
        A53=inhdr['A_5_3']
        A44=inhdr['A_4_4']
        A35=inhdr['A_3_5']
        A26=inhdr['A_2_6']
        A17=inhdr['A_1_7']
        A08=inhdr['A_0_8']

    if Aorder>=9:
        A90=inhdr['A_9_0']
        A81=inhdr['A_8_1']
        A72=inhdr['A_7_2']
        A63=inhdr['A_6_3']
        A54=inhdr['A_5_4']
        A45=inhdr['A_4_5']
        A36=inhdr['A_3_6']
        A27=inhdr['A_2_7']
        A18=inhdr['A_1_8']
        A09=inhdr['A_0_9']

    #coefficient B#
    if Border>=2:
        B20=inhdr['B_2_0']
        B11=inhdr['B_1_1']
        B02=inhdr['B_0_2']

    if Border>=3:
        B30=inhdr['B_3_0']
        B21=inhdr['B_2_1']
        B12=inhdr['B_1_2']
        B03=inhdr['B_0_3']

    if Border>=4:
        B40=inhdr['B_4_0']
        B31=inhdr['B_3_1']
        B22=inhdr['B_2_2']
        B13=inhdr['B_1_3']
        B04=inhdr['B_0_4']

    if Border>=5:
        B50=inhdr['B_5_0']
        B41=inhdr['B_4_1']
        B32=inhdr['B_3_2']
        B23=inhdr['B_2_3']
        B14=inhdr['B_1_4']
        B05=inhdr['B_0_5']

    if Border>=6:
        B60=inhdr['B_6_0']
        B51=inhdr['B_5_1']
        B42=inhdr['B_4_2']
        B33=inhdr['B_3_3']
        B24=inhdr['B_2_4']
        B15=inhdr['B_1_5']
        B06=inhdr['B_0_6']

    if Border>=7:
        B70=inhdr['B_7_0']
        B61=inhdr['B_6_1']
        B52=inhdr['B_5_2']
        B43=inhdr['B_4_3']
        B34=inhdr['B_3_4']
        B25=inhdr['B_2_5']
        B16=inhdr['B_1_6']
        B07=inhdr['B_0_7']

    if Border>=8:
        B80=inhdr['B_8_0']
        B71=inhdr['B_7_1']
        B62=inhdr['B_6_2']
        B53=inhdr['B_5_3']
        B44=inhdr['B_4_4']
        B35=inhdr['B_3_5']
        B26=inhdr['B_2_6']
        B17=inhdr['B_1_7']
        B08=inhdr['B_0_8']

    if Border>=9:
        B90=inhdr['B_9_0']
        B81=inhdr['B_8_1']
        B72=inhdr['B_7_2']
        B63=inhdr['B_6_3']
        B54=inhdr['B_5_4']
        B45=inhdr['B_4_5']
        B36=inhdr['B_3_6']
        B27=inhdr['B_2_7']
        B18=inhdr['B_1_8']
        B09=inhdr['B_0_9']

    # Now we compute the TPV outhdr (outhdr)
    #1st order#

    if Aorder>=1 or Border>=1:
        outhdr['PV1_0'] = 0.
        outhdr['PV1_1'] = 1.
        outhdr['PV1_2'] = 0.

        outhdr['PV2_0'] = 0.
        outhdr['PV2_1'] = 1.
        outhdr['PV2_2'] = 0.

    # r-terms
        outhdr['PV1_3'] = 0.
        outhdr['PV2_3'] = 0.

    #2nd order#

    if Aorder>=2 or Border>=2:

        D11=C11**2
        D12=C11*C21
        D13=C21**2

        D21=2*C11*C12
        D22=C11*C22+C12*C21
        D23=2*C21*C22

        D31=C12**2
        D32=C12*C22
        D33=C22**2

        D=np.matrix([[D11, D12, D13], [D21, D22, D23], [D31, D32, D33]], dtype=np.float64)
        ID=D.I

        XD1=A20*C11+B20*C12
        XD2=A11*C11+B11*C12
        XD3=A02*C11+B02*C12

        XD=np.matrix([[XD1], [XD2], [XD3]], dtype=np.float64)
        MD=np.dot(ID, XD)

        MD1=MD[0,0]
        MD2=MD[1,0]
        MD3=MD[2,0]

        outhdr['PV1_4'] = np.float64(MD1)
        outhdr['PV1_5'] = np.float64(MD2)
        outhdr['PV1_6'] = np.float64(MD3)

        YD1=A20*C21+B20*C22
        YD2=A11*C21+B11*C22
        YD3=A02*C21+B02*C22

        YD=np.matrix([[YD1], [YD2], [YD3]], dtype=np.float64)
        ND=np.dot(ID, YD)

        ND1=ND[2,0]
        ND2=ND[1,0]
        ND3=ND[0,0]

        outhdr['PV2_4'] = np.float64(ND1)
        outhdr['PV2_5'] = np.float64(ND2)
        outhdr['PV2_6'] = np.float64(ND3)

    #3rd order#

    if Aorder>=3 or Border>=3:

        E11=C11**3
        E12=(C11**2)*C21
        E13=C11*(C21**2)
        E14=C21**3

        E21=3*(C11**2)*C12
        E22=C11*(C11*C22+2*C12*C21)
        E23=C21*(2*C11*C22+C12*C21)
        E24=3*(C21**2)*C22

        E31=3*C11*(C12**2)
        E32=C12*(2*C11*C22+C12*C21)
        E33=C22*(C11*C22+2*C12*C21)
        E34=3*C21*(C22**2)

        E41=C12**3
        E42=(C12**2)*C22
        E43=C12*(C22**2)
        E44=C22**3

        E=np.matrix([[E11, E12, E13, E14], [E21, E22, E23, E24], [E31, E32, E33, E34], [E41, E42, E43, E44]], dtype=np.float64)
        IE=E.I

        XE1=A30*C11+B30*C12
        XE2=A21*C11+B21*C12
        XE3=A12*C11+B12*C12
        XE4=A03*C11+B03*C12

        XE=np.matrix([[XE1], [XE2], [XE3], [XE4]], dtype=np.float64)
        ME=np.dot(IE, XE)

        ME1=ME[0,0]
        ME2=ME[1,0]
        ME3=ME[2,0]
        ME4=ME[3,0]

        outhdr['PV1_7'] = np.float64(ME1)
        outhdr['PV1_8'] = np.float64(ME2)
        outhdr['PV1_9'] = np.float64(ME3)
        outhdr['PV1_10'] = np.float64(ME4)

        YE1=A30*C21+B30*C22
        YE2=A21*C21+B21*C22
        YE3=A12*C21+B12*C22
        YE4=A03*C21+B03*C22

        YE=np.matrix([[YE1], [YE2], [YE3], [YE4]], dtype=np.float64)
        NE=np.dot(IE, YE)

        NE1=NE[3,0]
        NE2=NE[2,0]
        NE3=NE[1,0]
        NE4=NE[0,0]

        outhdr['PV2_7'] = np.float64(NE1)
        outhdr['PV2_8'] = np.float64(NE2)
        outhdr['PV2_9'] = np.float64(NE3)
        outhdr['PV2_10'] = np.float64(NE4)


    # r^3-terms
        outhdr['PV1_11'] = 0.
        outhdr['PV2_11'] = 0.
        
    #4th order#

    if Aorder>=4 or Border>=4:

        F11=C11**4
        F12=(C11**3)*C21
        F13=(C11**2)*(C21**2)
        F14=C11*(C21**3)
        F15=C21**4

        F21=4*(C11**3)*C12
        F22=(C11**2)*(C11*C22+3*C12*C21)
        F23=2*C11*C21*(C11*C22+C12*C21)
        F24=(C21**2)*(3*C11*C22+C12*C21)
        F25=4*(C21**3)*C22

        F31=6*(C11**2)*(C12**2)
        F32=3*C11*C12*(C11*C22+C12*C21)
        F33=(C11**2)*(C22**2)+4*C11*C22*C12*C21+(C12**2)*(C21**2)
        F34=3*C21*C22*(C11*C22+C12*C21)
        F35=6*(C21**2)*(C22**2)

        F41=4*C11*(C12**3)
        F42=(C12**2)*(3*C11*C22+C12*C21)
        F43=2*C12*C22*(C11*C22+C12*C21)
        F44=(C22**2)*(C11*C22+3*C12*C21)
        F45=4*C21*(C22**3)

        F51=C12**4
        F52=(C12**3)*C22
        F53=(C12**2)*(C22**2)
        F54=C12*(C22**3)
        F55=C22**4

        F=np.matrix([[F11, F12, F13, F14, F15], [F21, F22, F23, F24, F25], [F31, F32, F33, F34, F35], [F41, F42, F43, F44, F45], [F51, F52, F53, F54, F55]], dtype=np.float64)
        IF=F.I

        XF1=A40*C11+B40*C12
        XF2=A31*C11+B31*C12
        XF3=A22*C11+B22*C12
        XF4=A13*C11+B13*C12
        XF5=A04*C11+B04*C12

        XF=np.matrix([[XF1], [XF2], [XF3], [XF4], [XF5]], dtype=np.float64)
        MF=np.dot(IF, XF)

        MF1=MF[0,0]
        MF2=MF[1,0]
        MF3=MF[2,0]
        MF4=MF[3,0]
        MF5=MF[4,0]

        outhdr['PV1_12'] = np.float64(MF1)
        outhdr['PV1_13'] = np.float64(MF2)
        outhdr['PV1_14'] = np.float64(MF3)
        outhdr['PV1_15'] = np.float64(MF4)
        outhdr['PV1_16'] = np.float64(MF5)

        YF1=A40*C21+B40*C22
        YF2=A31*C21+B31*C22
        YF3=A22*C21+B22*C22
        YF4=A13*C21+B13*C22
        YF5=A04*C21+B04*C22

        YF=np.matrix([[YF1], [YF2], [YF3], [YF4], [YF5]], dtype=np.float64)
        NF=np.dot(IF, YF)

        NF1=NF[4,0]
        NF2=NF[3,0]
        NF3=NF[2,0]
        NF4=NF[1,0]
        NF5=NF[0,0]

        outhdr['PV2_12'] = np.float64(NF1)
        outhdr['PV2_13'] = np.float64(NF2)
        outhdr['PV2_14'] = np.float64(NF3)
        outhdr['PV2_15'] = np.float64(NF4)
        outhdr['PV2_16'] = np.float64(NF5)


    #5th order#

    if Aorder>=5 or Border>=5:

        G11=C11**5
        G12=(C11**4)*C21
        G13=(C11**3)*(C21**2)
        G14=(C11**2)*(C21**3)
        G15=C11*(C21**4)
        G16=C21**5

        G21=5*(C11**4)*C12
        G22=(C11**3)*(C11*C22+4*C12*C21)
        G23=(C11**2)*C21*(2*C11*C22+3*C12*C21)
        G24=C11*(C21**2)*(3*C11*C22+2*C12*C21)
        G25=(C21**3)*(4*C11*C22+C12*C21)
        G26=5*(C21**4)*C22

        G31=10*(C11**3)*(C12**2)
        G32=2*(C11**2)*C12*(2*C11*C22+3*C12*C21)
        G33=C11*((C11**2)*(C22**2)+6*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        G34=C21*(3*(C11**2)*(C22**2)+6*C11*C22*C12*C21+(C12**2)*(C21**2))
        G35=2*(C21**2)*C22*(3*C11*C22+2*C12*C21)
        G36=10*(C21**3)*(C22**2)

        G41=10*(C11**2)*(C12**3)
        G42=2*C11*(C12**2)*(3*C11*C22+2*C12*C21)
        G43=C12*(3*(C11**2)*(C22**2)+6*C11*C22*C12*C21+(C12**2)*(C21**2))
        G44=C22*((C11**2)*(C22**2)+6*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        G45=2*C21*(C22**2)*(2*C11*C22+3*C12*C21)
        G46=10*(C21**2)*(C22**3)

        G51=5*C11*(C12**4)
        G52=(C12**3)*(4*C11*C22+C12*C21)
        G53=(C12**2)*C22*(3*C11*C22+2*C12*C21)
        G54=C12*(C22**2)*(2*C11*C22+3*C12*C21)
        G55=(C22**3)*(C11*C22+4*C12*C21)
        G56=5*C21*(C22**4)

        G61=C12**5
        G62=(C12**4)*C22
        G63=(C12**3)*(C22**2)
        G64=(C12**2)*(C22**3)
        G65=C12*(C22**4)
        G66=C22**5

        G=np.matrix([[G11, G12, G13, G14, G15, G16], [G21, G22, G23, G24, G25, G26], [G31, G32, G33, G34, G35, G36], [G41, G42, G43, G44, G45, G46], [G51, G52, G53, G54, G55, G56], [G61, G62, G63, G64, G65, G66]], dtype=np.float64)
        IG=G.I

        XG1=A50*C11+B50*C12
        XG2=A41*C11+B41*C12
        XG3=A32*C11+B32*C12
        XG4=A23*C11+B23*C12
        XG5=A14*C11+B14*C12
        XG6=A05*C11+B05*C12

        XG=np.matrix([[XG1], [XG2], [XG3], [XG4], [XG5], [XG6]], dtype=np.float64)
        MG=np.dot(IG, XG)

        MG1=MG[0,0]
        MG2=MG[1,0]
        MG3=MG[2,0]
        MG4=MG[3,0]
        MG5=MG[4,0]
        MG6=MG[5,0]

        outhdr['PV1_17'] = np.float64(MG1)
        outhdr['PV1_18'] = np.float64(MG2)
        outhdr['PV1_19'] = np.float64(MG3)
        outhdr['PV1_20'] = np.float64(MG4)
        outhdr['PV1_21'] = np.float64(MG5)
        outhdr['PV1_22'] = np.float64(MG6)

        YG1=A50*C21+B50*C22
        YG2=A41*C21+B41*C22
        YG3=A32*C21+B32*C22
        YG4=A23*C21+B23*C22
        YG5=A14*C21+B14*C22
        YG6=A05*C21+B05*C22

        YG=np.matrix([[YG1], [YG2], [YG3], [YG4], [YG5], [YG6]], dtype=np.float64)
        NG=np.dot(IG, YG)

        NG1=NG[5,0]
        NG2=NG[4,0]
        NG3=NG[3,0]
        NG4=NG[2,0]
        NG5=NG[1,0]
        NG6=NG[0,0]

        outhdr['PV2_17'] = np.float64(NG1)
        outhdr['PV2_18'] = np.float64(NG2)
        outhdr['PV2_19'] = np.float64(NG3)
        outhdr['PV2_20'] = np.float64(NG4)
        outhdr['PV2_21'] = np.float64(NG5)
        outhdr['PV2_22'] = np.float64(NG6)

    # r^5-terms
        outhdr['PV1_23'] = 0.
        outhdr['PV2_23'] = 0.
            
    #6th order#

    if Aorder>=6 or Border>=6:

        H11=C11**6
        H12=(C11**5)*C21
        H13=(C11**4)*(C21**2)
        H14=(C11**3)*(C21**3)
        H15=(C11**2)*(C21**4)
        H16=C11*(C21**5)
        H17=C21**6

        H21=6*(C11**5)*C12
        H22=(C11**4)*(C11*C22+5*C12*C21)
        H23=(C11**3)*C21*(2*C11*C22+4*C12*C21)
        H24=3*(C11**2)*(C21**2)*(C11*C22+C12*C21)
        H25=C11*(C21**3)*(4*C11*C22+2*C12*C21)
        H26=(C21**4)*(5*C11*C22+C12*C21)
        H27=6*(C21**5)*C22

        H31=15*(C11**4)*(C12**2)
        H32=5*(C11**3)*C12*(C11*C22+2*C12*C21)
        H33=(C11**2)*((C11**2)*(C22**2)+8*C11*C22*C12*C21+6*(C12**2)*(C21**2))
        H34=3*C11*C21*((C11**2)*(C22**2)+3*C11*C22*C12*C21+(C12**2)*(C21**2))
        H35=(C21**2)*(6*(C11**2)*(C22**2)+8*C11*C22*C12*C21+(C12**2)*(C21**2))
        H36=5*(C21**3)*C22*(2*C11*C22+C12*C21)
        H37=15*(C21**4)*(C22**2)

        H41=20*(C11**3)*(C12**3)
        H42=10*(C11**2)*(C12**2)*(C11*C22+C12*C21)
        H43=4*C11*C12*((C11**2)*(C22**2)+3*C11*C22*C12*C21+(C12**2)*(C21**2))
        H44=(C11**3)*(C22**3)+9*(C11**2)*(C22**2)*C12*C21+9*C11*C22*(C12**2)*(C21**2)+(C12**3)*(C21**3)
        H45=4*C21*C22*((C11**2)*(C22**2)+3*C11*C22*C12*C21+(C12**2)*(C21**2))
        H46=10*(C21**2)*(C22**2)*(C11*C22+C12*C21)
        H47=20*(C21**3)*(C22**3)

        H51=15*(C11**2)*(C12**4)
        H52=5*C11*(C12**3)*(2*C11*C22+C12*C21)
        H53=(C12**2)*(6*(C11**2)*(C22**2)+8*C11*C22*C12*C21+(C12**2)*(C21**2))
        H54=3*C12*C22*((C11**2)*(C22**2)+3*C11*C22*C12*C21+(C12**2)*(C21**2))
        H55=(C22**2)*((C11**2)*(C22**2)+8*C11*C22*C12*C21+6*(C12**2)*(C21**2))
        H56=5*C21*(C22**3)*(C11*C22+2*C12*C21)
        H57=15*(C21**2)*(C22**4)

        H61=6*C11*(C12**5)
        H62=(C12**4)*(5*C11*C22+C12*C21)
        H63=(C12**3)*C22*(4*C11*C22+2*C12*C21)
        H64=3*(C12**2)*(C22**2)*(C11*C22+C12*C21)
        H65=C12*(C22**3)*(2*C11*C22+4*C12*C21)
        H66=(C22**4)*(C11*C22+5*C12*C21)
        H67=6*C21*(C22**5)

        H71=C12**6
        H72=(C12**5)*C22
        H73=(C12**4)*(C22**2)
        H74=(C12**3)*(C22**3)
        H75=(C12**2)*(C22**4)
        H76=C12*(C22**5)
        H77=C22**6

        H=np.matrix([[H11, H12, H13, H14, H15, H16, H17], [H21, H22, H23, H24, H25, H26, H27], [H31, H32, H33, H34, H35, H36, H37], [H41, H42, H43, H44, H45, H46, H47], [H51, H52, H53, H54, H55, H56, H57], [H61, H62, H63, H64, H65, H66, H67], [H71, H72, H73, H74, H75, H76, H77]], dtype=np.float64)
        IH=H.I

        XH1=A60*C11+B60*C12
        XH2=A51*C11+B51*C12
        XH3=A42*C11+B42*C12
        XH4=A33*C11+B33*C12
        XH5=A24*C11+B24*C12
        XH6=A15*C11+B15*C12
        XH7=A06*C11+B06*C12

        XH=np.matrix([[XH1], [XH2], [XH3], [XH4], [XH5], [XH6], [XH7]], dtype=np.float64)
        MH=np.dot(IH, XH)

        MH1=MH[0,0]
        MH2=MH[1,0]
        MH3=MH[2,0]
        MH4=MH[3,0]
        MH5=MH[4,0]
        MH6=MH[5,0]
        MH7=MH[6,0]

        outhdr['PV1_24'] = np.float64(MH1)
        outhdr['PV1_25'] = np.float64(MH2)
        outhdr['PV1_26'] = np.float64(MH3)
        outhdr['PV1_27'] = np.float64(MH4)
        outhdr['PV1_28'] = np.float64(MH5)
        outhdr['PV1_29'] = np.float64(MH6)
        outhdr['PV1_30'] = np.float64(MH7)

        YH1=A60*C21+B60*C22
        YH2=A51*C21+B51*C22
        YH3=A42*C21+B42*C22
        YH4=A33*C21+B33*C22
        YH5=A24*C21+B24*C22
        YH6=A15*C21+B15*C22
        YH7=A06*C21+B06*C22

        YH=np.matrix([[YH1], [YH2], [YH3], [YH4], [YH5], [YH6], [YH7]], dtype=np.float64)
        NH=np.dot(IH, YH)

        NH1=NH[6,0]
        NH2=NH[5,0]
        NH3=NH[4,0]
        NH4=NH[3,0]
        NH5=NH[2,0]
        NH6=NH[1,0]
        NH7=NH[0,0]

        outhdr['PV2_24'] = np.float64(NH1)
        outhdr['PV2_25'] = np.float64(NH2)
        outhdr['PV2_26'] = np.float64(NH3)
        outhdr['PV2_27'] = np.float64(NH4)
        outhdr['PV2_28'] = np.float64(NH5)
        outhdr['PV2_29'] = np.float64(NH6)
        outhdr['PV2_30'] = np.float64(NH7)


    #7th order#

    if Aorder>=7 or Border>=7:

        J11=C11**7
        J12=(C11**6)*C21
        J13=(C11**5)*(C21**2)
        J14=(C11**4)*(C21**3)
        J15=(C11**3)*(C21**4)
        J16=(C11**2)*(C21**5)
        J17=C11*(C21**6)
        J18=C21**7

        J21=7*(C11**6)*C12
        J22=(C11**5)*(C11*C22+6*C12*C21)
        J23=(C11**4)*C21*(2*C11*C22+5*C12*C21)
        J24=(C11**3)*(C21**2)*(3*C11*C22+4*C12*C21)
        J25=(C11**2)*(C21**3)*(4*C11*C22+3*C12*C21)
        J26=C11*(C21**4)*(5*C11*C22+2*C12*C21)
        J27=(C21**5)*(6*C11*C22+C12*C21)
        J28=7*(C21**6)*C22

        J31=21*(C11**5)*(C12**2)
        J32=3*(C11**4)*C12*(2*C11*C22+5*C12*C21)
        J33=(C11**3)*((C11**2)*(C22**2)+10*C11*C22*C12*C21+10*(C12**2)*(C21**2))
        J34=3*(C11**2)*C21*((C11**2)*(C22**2)+4*C11*C22*C12*C21+2*(C12**2)*(C21**2))
        J35=3*C11*(C21**2)*(2*(C11**2)*(C22**2)+4*C11*C22*C12*C21+(C12**2)*(C21**2))
        J36=(C21**3)*(10*(C11**2)*(C22**2)+10*C11*C22*C12*C21+(C12**2)*(C21**2))
        J37=3*(C21**4)*C22*(5*C11*C22+2*C12*C21)
        J38=21*(C21**5)*(C22**2)

        J41=35*(C11**4)*(C12**3)
        J42=5*(C11**3)*(C12**2)*(3*C11*C22+4*C12*C21)
        J43=5*(C11**2)*C12*((C11**2)*(C22**2)+4*C11*C22*C12*C21+2*(C12**2)*(C21**2))
        J44=C11*((C11**3)*(C22**3)+12*(C11**2)*(C22**2)*C12*C21+18*C11*C22*(C12**2)*(C21**2)+4*(C12**3)*(C21**3))
        J45=C21*(4*(C11**3)*(C22**3)+18*(C11**2)*(C22**2)*C12*C21+12*C11*C22*(C12**2)*(C21**2)+(C12**3)*(C21**3))
        J46=5*(C21**2)*C22*(2*(C11**2)*(C22**2)+4*C11*C22*C12*C21+(C12**2)*(C21**2))
        J47=5*(C21**3)*(C22**2)*(4*C11*C22+3*C12*C21)
        J48=35*(C21**4)*(C22**3)

        J51=35*(C11**3)*(C12**4)
        J52=5*(C11**2)*(C12**3)*(4*C11*C22+3*C12*C21)
        J53=5*C11*(C12**2)*(2*(C11**2)*(C22**2)+4*C11*C22*C12*C21+(C12**2)*(C21**2))
        J54=C12*(4*(C11**3)*(C22**3)+18*(C11**2)*(C22**2)*C12*C21+12*C11*C22*(C12**2)*(C21**2)+(C12**3)*(C21**3))
        J55=C22*((C11**3)*(C22**3)+12*(C11**2)*(C22**2)*C12*C21+18*C11*C22*(C12**2)*(C21**2)+4*(C12**3)*(C21**3))
        J56=5*C21*(C22**2)*((C11**2)*(C22**2)+4*C11*C22*C12*C21+2*(C12**2)*(C21**2))
        J57=5*(C21**2)*(C22**3)*(3*C11*C22+4*C12*C21)
        J58=35*(C21**3)*(C22**4)

        J61=21*(C11**2)*(C12**5)
        J62=3*C11*(C12**4)*(5*C11*C22+2*C12*C21)
        J63=(C12**3)*(10*(C11**2)*(C22**2)+10*C11*C22*C12*C21+(C12**2)*(C21**2))
        J64=3*(C12**2)*C22*(2*(C11**2)*(C22**2)+4*C11*C22*C12*C21+(C12**2)*(C21**2))
        J65=3*C12*(C22**2)*((C11**2)*(C22**2)+4*C11*C22*C12*C21+2*(C12**2)*(C21**2))
        J66=(C22**3)*((C11**2)*(C22**2)+10*C11*C22*C12*C21+10*(C12**2)*(C21**2))
        J67=3*C21*(C22**4)*(2*C11*C22+5*C12*C21)
        J68=21*(C21**2)*(C22**5)

        J71=7*C11*(C12**6)
        J72=(C12**5)*(6*C11*C22+C12*C21)
        J73=(C12**4)*C22*(5*C11*C22+2*C12*C21)
        J74=(C12**3)*(C22**2)*(4*C11*C22+3*C12*C21)
        J75=(C12**2)*(C22**3)*(3*C11*C22+4*C12*C21)
        J76=C12*(C22**4)*(2*C11*C22+5*C12*C21)
        J77=(C22**5)*(C11*C22+6*C12*C21)
        J78=7*C21*(C22**6)

        J81=C12**7
        J82=(C12**6)*C22
        J83=(C12**5)*(C22**2)
        J84=(C12**4)*(C22**3)
        J85=(C12**3)*(C22**4)
        J86=(C12**2)*(C22**5)
        J87=C12*(C22**6)
        J88=C22**7

        J=np.matrix([[J11, J12, J13, J14, J15, J16, J17, J18], [J21, J22, J23, J24, J25, J26, J27, J28], [J31, J32, J33, J34, J35, J36, J37, J38], [J41, J42, J43, J44, J45, J46, J47, J48], [J51, J52, J53, J54, J55, J56, J57, J58], [J61, J62, J63, J64, J65, J66, J67, J68], [J71, J72, J73, J74, J75, J76, J77, J78], [J81, J82, J83, J84, J85, J86, J87, J88]], dtype=np.float64)
        IJ=J.I

        XJ1=A70*C11+B70*C12
        XJ2=A61*C11+B61*C12
        XJ3=A52*C11+B52*C12
        XJ4=A43*C11+B43*C12
        XJ5=A34*C11+B34*C12
        XJ6=A25*C11+B25*C12
        XJ7=A16*C11+B16*C12
        XJ8=A07*C11+B07*C12

        XJ=np.matrix([[XJ1], [XJ2], [XJ3], [XJ4], [XJ5], [XJ6], [XJ7], [XJ8]], dtype=np.float64)
        MJ=np.dot(IJ, XJ)

        MJ1=MJ[0,0]
        MJ2=MJ[1,0]
        MJ3=MJ[2,0]
        MJ4=MJ[3,0]
        MJ5=MJ[4,0]
        MJ6=MJ[5,0]
        MJ7=MJ[6,0]
        MJ8=MJ[7,0]

        outhdr['PV1_31'] = np.float64(MJ1)
        outhdr['PV1_32'] = np.float64(MJ2)
        outhdr['PV1_33'] = np.float64(MJ3)
        outhdr['PV1_34'] = np.float64(MJ4)
        outhdr['PV1_35'] = np.float64(MJ5)
        outhdr['PV1_36'] = np.float64(MJ6)
        outhdr['PV1_37'] = np.float64(MJ7)
        outhdr['PV1_38'] = np.float64(MJ8)


        YJ1=A70*C21+B70*C22
        YJ2=A61*C21+B61*C22
        YJ3=A52*C21+B52*C22
        YJ4=A43*C21+B43*C22
        YJ5=A34*C21+B34*C22
        YJ6=A25*C21+B25*C22
        YJ7=A16*C21+B16*C22
        YJ8=A07*C21+B07*C22

        YJ=np.matrix([[YJ1], [YJ2], [YJ3], [YJ4], [YJ5], [YJ6], [YJ7], [YJ8]], dtype=np.float64)
        NJ=np.dot(IJ, YJ)

        NJ1=NJ[7,0]
        NJ2=NJ[6,0]
        NJ3=NJ[5,0]
        NJ4=NJ[4,0]
        NJ5=NJ[3,0]
        NJ6=NJ[2,0]
        NJ7=NJ[1,0]
        NJ8=NJ[0,0]

        outhdr['PV2_31'] = np.float64(NJ1)
        outhdr['PV2_32'] = np.float64(NJ2)
        outhdr['PV2_33'] = np.float64(NJ3)
        outhdr['PV2_34'] = np.float64(NJ4)
        outhdr['PV2_35'] = np.float64(NJ5)
        outhdr['PV2_36'] = np.float64(NJ6)
        outhdr['PV2_37'] = np.float64(NJ7)
        outhdr['PV2_38'] = np.float64(NJ8)

    # r^7-terms
        outhdr['PV1_39'] = 0.
        outhdr['PV2_39'] = 0.

    #8th order#

    if Aorder>=8 or Border>=8:

        P11=C11**8
        P12=(C11**7)*C21
        P13=(C11**6)*(C21**2)
        P14=(C11**5)*(C21**3)
        P15=(C11**4)*(C21**4)
        P16=(C11**3)*(C21**5)
        P17=(C11**2)*(C21**6)
        P18=C11*(C21**7)
        P19=C21**8

        P21=8*(C11**7)*C12
        P22=(C11**6)*(C11*C22+7*C12*C21)
        P23=2*(C11**5)*C21*(C11*C22+3*C12*C21)
        P24=(C11**4)*(C21**2)*(3*C11*C22+5*C12*C21)
        P25=4*(C11**3)*(C21**3)*(C11*C22+C12*C21)
        P26=(C11**2)*(C21**4)*(5*C11*C22+3*C12*C21)
        P27=2*C11*(C21**5)*(3*C11*C22+C12*C21)
        P28=(C21**6)*(7*C11*C22+C12*C21)
        P29=8*(C21**7)*C22

        P31=28*(C11**6)*(C12**2)
        P32=7*(C11**5)*C12*(C11*C22+3*C12*C21)
        P33=(C11**4)*((C11**2)*(C22**2)+12*C11*C22*C12*C21+15*(C12**2)*(C21**2))
        P34=(C11**3)*C21*(3*(C11**2)*(C22**2)+15*C11*C22*C12*C21+10*(C12**2)*(C21**2))
        P35=2*(C11**2)*(C21**2)*(3*(C11**2)*(C22**2)+8*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        P36=C11*(C21**3)*(10*(C11**2)*(C22**2)+15*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        P37=(C21**4)*(15*(C11**2)*(C22**2)+12*C11*C22*C12*C21+(C12**2)*(C21**2))
        P38=7*(C21**5)*C22*(3*C11*C22+C12*C21)
        P39=28*(C21**6)*(C22**2)

        P41=56*(C11**5)*(C12**3)
        P42=7*(C11**4)*(C12**2)*(3*C11*C22+5*C12*C21)
        P43=2*(C11**3)*C12*(3*(C11**2)*(C22**2)+15*C11*C22*C12*C21+10*(C12**2)*(C21**2))
        P44=(C11**2)*((C11**3)*(C22**3)+15*(C11**2)*(C22**2)*C12*C21+30*C11*C22*(C12**2)*(C21**2)+10*(C12**3)*(C21**3))
        P45=4*C11*C21*((C11**3)*(C22**3)+6*(C11**2)*(C22**2)*C12*C21+6*C11*C22*(C12**2)*(C21**2)+(C12**3)*(C21**3))
        P46=(C21**2)*(10*(C11**3)*(C22**3)+30*(C11**2)*(C22**2)*C12*C21+15*C11*C22*(C12**2)*(C21**2)+(C12**3)*(C21**3))
        P47=2*(C21**3)*C22*(10*(C11**2)*(C22**2)+15*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        P48=7*(C21**4)*(C22**2)*(5*C11*C22+3*C12*C21)
        P49=56*(C21**5)*(C22**3)

        P51=70*(C11**4)*(C12**4)
        P52=35*(C11**3)*(C12**3)*(C11*C22+C12*C21)
        P53=5*(C11**2)*(C12**2)*(3*(C11**2)*(C22**2)+8*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        P54=5*C11*C12*((C11**3)*(C22**3)+6*(C11**2)*(C22**2)*C12*C21+6*C11*C22*(C12**2)*(C21**2)+(C12**3)*(C21**3))
        P55=(C11**4)*(C22**4)+16*(C11**3)*(C22**3)*C12*C21+36*(C11**2)*(C22**2)*(C12**2)*(C21**2)+16*C11*C22*(C12**3)*(C21**3)+(C12**4)*(C21**4)
        P56=5*C21*C22*((C11**3)*(C22**3)+6*(C11**2)*(C22**2)*C12*C21+6*C11*C22*(C12**2)*(C21**2)+(C12**3)*(C21**3))
        P57=5*(C21**2)*(C22**2)*(3*(C11**2)*(C22**2)+8*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        P58=35*(C21**3)*(C22**3)*(C11*C22+C12*C21)
        P59=70*(C21**4)*(C22**4)

        P61=56*(C11**3)*(C12**5)
        P62=7*(C11**2)*(C12**4)*(5*C11*C22+3*C12*C21)
        P63=2*C11*(C12**3)*(10*(C11**2)*(C22**2)+15*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        P64=(C12**2)*(10*(C11**3)*(C22**3)+30*(C11**2)*(C22**2)*C12*C21+15*C11*C22*(C12**2)*(C21**2)+(C12**3)*(C21**3))
        P65=4*C12*C22*((C11**3)*(C22**3)+6*(C11**2)*(C22**2)*C12*C21+6*C11*C22*(C12**2)*(C21**2)+(C12**3)*(C21**3))
        P66=(C22**2)*((C11**3)*(C22**3)+15*(C11**2)*(C22**2)*C12*C21+30*C11*C22*(C12**2)*(C21**2)+10*(C12**3)*(C21**3))
        P67=2*C21*(C22**3)*(3*(C11**2)*(C22**2)+15*C11*C22*C12*C21+10*(C12**2)*(C21**2))
        P68=7*(C21**2)*(C22**4)*(3*C11*C22+5*C12*C21)
        P69=56*(C21**3)*(C22**5)

        P71=28*(C11**2)*(C12**6)
        P72=7*C11*(C12**5)*(3*C11*C22+C12*C21)
        P73=(C12**4)*(15*(C11**2)*(C22**2)+12*C11*C22*C12*C21+(C12**2)*(C21**2))
        P74=(C12**3)*C22*(10*(C11**2)*(C22**2)+15*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        P75=2*(C12**2)*(C22**2)*(3*(C11**2)*(C22**2)+8*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        P76=C12*(C22**3)*(3*(C11**2)*(C22**2)+15*C11*C22*C12*C21+10*(C12**2)*(C21**2))
        P77=(C22**4)*((C11**2)*(C22**2)+12*C11*C22*C12*C21+15*(C12**2)*(C21**2))
        P78=7*C21*(C22**5)*(C11*C22+3*C12*C21)
        P79=28*(C21**2)*(C22**6)

        P81=8*C11*(C12**7)
        P82=(C12**6)*(7*C11*C22+C12*C21)
        P83=2*(C12**5)*C22*(3*C11*C22+C12*C21)
        P84=(C12**4)*(C22**2)*(5*C11*C22+3*C12*C21)
        P85=4*(C12**3)*(C22**3)*(C11*C22+C12*C21)
        P86=(C12**2)*(C22**4)*(3*C11*C22+5*C12*C21)
        P87=2*C12*(C22**5)*(C11*C22+3*C12*C21)
        P88=(C22**6)*(C11*C22+7*C12*C21)
        P89=8*C21*(C22**7)

        P91=C12**8
        P92=(C12**7)*C22
        P93=(C12**6)*(C22**2)
        P94=(C12**5)*(C22**3)
        P95=(C12**4)*(C22**4)
        P96=(C12**3)*(C22**5)
        P97=(C12**2)*(C22**6)
        P98=C12*(C22**7)
        P99=C22**8

        P=np.matrix([[P11, P12, P13, P14, P15, P16, P17, P18, P19], [P21, P22, P23, P24, P25, P26, P27, P28, P29], [P31, P32, P33, P34, P35, P36, P37, P38, P39], [P41, P42, P43, P44, P45, P46, P47, P48, P49], [P51, P52, P53, P54, P55, P56, P57, P58, P59], [P61, P62, P63, P64, P65, P66, P67, P68, P69], [P71, P72, P73, P74, P75, P76, P77, P78, P79], [P81, P82, P83, P84, P85, P86, P87, P88, P89], [P91, P92, P93, P94, P95, P96, P97, P98, P99]], dtype=np.float)
        IP=P.I

        XP1=A80*C11+B80*C12
        XP2=A71*C11+B71*C12
        XP3=A62*C11+B62*C12
        XP4=A53*C11+B53*C12
        XP5=A44*C11+B44*C12
        XP6=A35*C11+B35*C12
        XP7=A26*C11+B26*C12
        XP8=A17*C11+B17*C12
        XP9=A08*C11+B08*C12

        XP=np.matrix([[XP1], [XP2], [XP3], [XP4], [XP5], [XP6], [XP7], [XP8], [XP9]], dtype=np.float)
        MP=np.dot(IP, XP)

        MP1=MP[0,0]
        MP2=MP[1,0]
        MP3=MP[2,0]
        MP4=MP[3,0]
        MP5=MP[4,0]
        MP6=MP[5,0]
        MP7=MP[6,0]
        MP8=MP[7,0]
        MP9=MP[8,0]

        outhdr['PV1_40'] = np.float64(MP1)
        outhdr['PV1_41'] = np.float64(MP2)
        outhdr['PV1_42'] = np.float64(MP3)
        outhdr['PV1_43'] = np.float64(MP4)
        outhdr['PV1_44'] = np.float64(MP5)
        outhdr['PV1_45'] = np.float64(MP6)
        outhdr['PV1_46'] = np.float64(MP7)
        outhdr['PV1_47'] = np.float64(MP8)
        outhdr['PV1_48'] = np.float64(MP9)


        YP1=A80*C21+B80*C22
        YP2=A71*C21+B71*C22
        YP3=A62*C21+B62*C22
        YP4=A53*C21+B53*C22
        YP5=A44*C21+B44*C22
        YP6=A35*C21+B35*C22
        YP7=A26*C21+B26*C22
        YP8=A17*C21+B17*C22
        YP9=A08*C21+B08*C22

        YP=np.matrix([[YP1], [YP2], [YP3], [YP4], [YP5], [YP6], [YP7], [YP8], [YP9]], dtype=np.float)
        NP=np.dot(IP, YP)

        NP1=NP[8,0]
        NP2=NP[7,0]
        NP3=NP[6,0]
        NP4=NP[5,0]
        NP5=NP[4,0]
        NP6=NP[3,0]
        NP7=NP[2,0]
        NP8=NP[1,0]
        NP9=NP[0,0]

        outhdr['PV2_40'] = np.float64(NP1)
        outhdr['PV2_41'] = np.float64(NP2)
        outhdr['PV2_42'] = np.float64(NP3)
        outhdr['PV2_43'] = np.float64(NP4)
        outhdr['PV2_44'] = np.float64(NP5)
        outhdr['PV2_45'] = np.float64(NP6)
        outhdr['PV2_46'] = np.float64(NP7)
        outhdr['PV2_47'] = np.float64(NP8)
        outhdr['PV2_48'] = np.float64(NP9)

    #9th order#

    if Aorder>=9 or Border>=9:

        Q11=C11**9
        Q12=(C11**8)*C21
        Q13=(C11**7)*(C21**2)
        Q14=(C11**6)*(C21**3)
        Q15=(C11**5)*(C21**4)
        Q16=(C11**4)*(C21**5)
        Q17=(C11**3)*(C21**6)
        Q18=(C11**2)*(C21**7)
        Q19=C11*(C21**8)
        Q10=C21**9

        Q21=9*(C11**8)*C12
        Q22=(C11**7)*(C11*C22+8*C12*C21)
        Q23=(C11**6)*C21*(2*C11*C22+7*C12*C21)
        Q24=3*(C11**5)*(C21**2)*(C11*C22+2*C12*C21)
        Q25=(C11**4)*(C21**3)*(4*C11*C22+5*C12*C21)
        Q26=(C11**3)*(C21**4)*(5*C11*C22+4*C12*C21)
        Q27=3*(C11**2)*(C21**5)*(2*C11*C22+C12*C21)
        Q28=C11*(C21**6)*(7*C11*C22+2*C12*C21)
        Q29=(C21**7)*(8*C11*C22+C12*C21)
        Q20=9*(C21**8)*C22

        Q31=36*(C11**7)*(C12**2)
        Q32=4*(C11**6)*C12*(2*C11*C22+7*C12*C21)
        Q33=(C11**5)*((C11**2)*(C22**2)+14*C11*C22*C12*C21+21*(C12**2)*(C21**2))
        Q34=3*(C11**4)*C21*((C11**2)*(C22**2)+6*C11*C22*C12*C21+5*(C12**2)*(C21**2))
        Q35=2*(C11**3)*(C21**2)*(3*(C11**2)*(C22**2)+10*C11*C22*C12*C21+5*(C12**2)*(C21**2))
        Q36=2*(C11**2)*(C21**3)*(5*(C11**2)*(C22**2)+10*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        Q37=3*C11*(C21**4)*(5*(C11**2)*(C22**2)+6*C11*C22*C12*C21+(C12**2)*(C21**2))
        Q38=(C21**5)*(21*(C11**2)*(C22**2)+14*C11*C22*C12*C21+(C12**2)*(C21**2))
        Q39=4*(C21**6)*C22*(7*C11*C22+2*C12*C21)
        Q30=36*(C21**7)*(C22**2)

        Q41=84*(C11**6)*(C12**3)
        Q42=28*(C11**5)*(C12**2)*(C11*C22+2*C12*C21)
        Q43=7*(C11**4)*C12*((C11**2)*(C22**2)+6*C11*C22*C12*C21+5*(C12**2)*(C21**2))
        Q44=(C11**3)*((C11**3)*(C22**3)+18*(C11**2)*(C22**2)*C12*C21+45*C11*C22*(C12**2)*(C21**2)+20*(C12**3)*(C21**3))
        Q45=2*(C11**2)*C21*(2*(C11**3)*(C22**3)+15*(C11**2)*(C22**2)*C12*C21+20*C11*C22*(C12**2)*(C21**2)+5*(C12**3)*(C21**3))
        Q46=2*C11*(C21**2)*(5*(C11**3)*(C22**3)+20*(C11**2)*(C22**2)*C12*C21+15*C11*C22*(C12**2)*(C21**2)+2*(C12**3)*(C21**3))
        Q47=(C21**3)*(20*(C11**3)*(C22**3)+45*(C11**2)*(C22**2)*C12*C21+18*C11*C22*(C12**2)*(C21**2)+(C12**3)*(C21**3))
        Q48=7*(C21**4)*C22*(5*(C11**2)*(C22**2)+6*C11*C22*C12*C21+(C12**2)*(C21**2))
        Q49=28*(C21**5)*(C22**2)*(2*C11*C22+C12*C21)
        Q40=84*(C21**6)*(C22**3)

        Q51=126*(C11**5)*(C12**4)
        Q52=14*(C11**4)*(C12**3)*(4*C11*C22+5*C12*C21)
        Q53=7*(C11**3)*(C12**2)*(3*(C11**2)*(C22**2)+10*C11*C22*C12*C21+5*(C12**2)*(C21**2))
        Q54=3*(C11**2)*C12*(2*(C11**3)*(C22**3)+15*(C11**2)*(C22**2)*C12*C21+20*C11*C22*(C12**2)*(C21**2)+5*(C12**3)*(C21**3))
        Q55=C11*((C11**4)*(C22**4)+20*(C11**3)*(C22**3)*C12*C21+60*(C11**2)*(C22**2)*(C12**2)*(C21**2)+40*C11*C22*(C12**3)*(C21**3)+5*(C12**4)*(C21**4))
        Q56=C21*(5*(C11**4)*(C22**4)+40*(C11**3)*(C22**3)*C12*C21+60*(C11**2)*(C22**2)*(C12**2)*(C21**2)+20*C11*C22*(C12**3)*(C21**3)+(C12**4)*(C21**4))
        Q57=3*(C21**2)*C22*(5*(C11**3)*(C22**3)+20*(C11**2)*(C22**2)*C12*C21+15*C11*C22*(C12**2)*(C21**2)+2*(C12**3)*(C21**3))
        Q58=7*(C21**3)*(C22**2)*(5*(C11**2)*(C22**2)+10*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        Q59=14*(C21**4)*(C22**3)*(5*C11*C22+4*C12*C21)
        Q50=126*(C21**5)*(C22**4)

        Q61=126*(C11**4)*(C12**5)
        Q62=14*(C11**3)*(C12**4)*(5*C11*C22+4*C12*C21)
        Q63=7*(C11**2)*(C12**3)*(5*(C11**2)*(C22**2)+10*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        Q64=3*C11*(C12**2)*(5*(C11**3)*(C22**3)+20*(C11**2)*(C22**2)*C12*C21+15*C11*C22*(C12**2)*(C21**2)+2*(C12**3)*(C21**3))
        Q65=C12*(5*(C11**4)*(C22**4)+40*(C11**3)*(C22**3)*C12*C21+60*(C11**2)*(C22**2)*(C12**2)*(C21**2)+20*C11*C22*(C12**3)*(C21**3)+(C12**4)*(C21**4))
        Q66=C22*((C11**4)*(C22**4)+20*(C11**3)*(C22**3)*C12*C21+60*(C11**2)*(C22**2)*(C12**2)*(C21**2)+40*C11*C22*(C12**3)*(C21**3)+5*(C12**4)*(C21**4))
        Q67=3*C21*(C22**2)*(2*(C11**3)*(C22**3)+15*(C11**2)*(C22**2)*C12*C21+20*C11*C22*(C12**2)*(C21**2)+5*(C12**3)*(C21**3))
        Q68=7*(C21**2)*(C22**3)*(3*(C11**2)*(C22**2)+10*C11*C22*C12*C21+5*(C12**2)*(C21**2))
        Q69=14*(C21**3)*(C22**4)*(4*C11*C22+5*C12*C21)
        Q60=126*(C21**4)*(C22**5)

        Q71=84*(C11**3)*(C12**6)
        Q72=28*(C11**2)*(C12**5)*(2*C11*C22+C12*C21)
        Q73=7*C11*(C12**4)*(5*(C11**2)*(C22**2)+6*C11*C22*C12*C21+(C12**2)*(C21**2))
        Q74=(C12**3)*(20*(C11**3)*(C22**3)+45*(C11**2)*(C22**2)*C12*C21+18*C11*C22*(C12**2)*(C21**2)+(C12**3)*(C21**3))
        Q75=2*(C12**2)*C22*(5*(C11**3)*(C22**3)+20*(C11**2)*(C22**2)*C12*C21+15*C11*C22*(C12**2)*(C21**2)+2*(C12**3)*(C21**3))
        Q76=2*C12*(C22**2)*(2*(C11**3)*(C22**3)+15*(C11**2)*(C22**2)*C12*C21+20*C11*C22*(C12**2)*(C21**2)+5*(C12**3)*(C21**3))
        Q77=(C22**3)*((C11**3)*(C22**3)+18*(C11**2)*(C22**2)*C12*C21+45*C11*C22*(C12**2)*(C21**2)+20*(C12**3)*(C21**3))
        Q78=7*C21*(C22**4)*((C11**2)*(C22**2)+6*C11*C22*C12*C21+5*(C12**2)*(C21**2))
        Q79=28*(C21**2)*(C22**5)*(C11*C22+2*C12*C21)
        Q70=84*(C21**3)*(C22**6)

        Q81=36*(C11**2)*(C12**7)
        Q82=4*C11*(C12**6)*(7*C11*C22+2*C12*C21)
        Q83=(C12**5)*(21*(C11**2)*(C22**2)+14*C11*C22*C12*C21+(C12**2)*(C21**2))
        Q84=3*(C12**4)*C22*(5*(C11**2)*(C22**2)+6*C11*C22*C12*C21+(C12**2)*(C21**2))
        Q85=2*(C12**3)*(C22**2)*(5*(C11**2)*(C22**2)+10*C11*C22*C12*C21+3*(C12**2)*(C21**2))
        Q86=2*(C12**2)*(C22**3)*(3*(C11**2)*(C22**2)+10*C11*C22*C12*C21+5*(C12**2)*(C21**2))
        Q87=3*C12*(C22**4)*((C11**2)*(C22**2)+6*C11*C22*C12*C21+5*(C12**2)*(C21**2))
        Q88=(C22**5)*((C11**2)*(C22**2)+14*C11*C22*C12*C21+21*(C12**2)*(C21**2))
        Q89=4*C21*(C22**6)*(2*C11*C22+7*C12*C21)
        Q80=36*(C21**2)*(C22**7)

        Q91=9*C11*(C12**8)
        Q92=(C12**7)*(8*C11*C22+C12*C21)
        Q93=(C12**6)*C22*(7*C11*C22+2*C12*C21)
        Q94=3*(C12**5)*(C22**2)*(2*C11*C22+C12*C21)
        Q95=(C12**4)*(C22**3)*(5*C11*C22+4*C12*C21)
        Q96=(C12**3)*(C22**4)*(4*C11*C22+5*C12*C21)
        Q97=3*(C12**2)*(C22**5)*(C11*C22+2*C12*C21)
        Q98=C12*(C22**6)*(2*C11*C22+7*C12*C21)
        Q99=(C22**7)*(C11*C22+8*C12*C21)
        Q90=9*C21*(C22**8)

        Q01=C12**9
        Q02=(C12**8)*C22
        Q03=(C12**7)*(C22**2)
        Q04=(C12**6)*(C22**3)
        Q05=(C12**5)*(C22**4)
        Q06=(C12**4)*(C22**5)
        Q07=(C12**3)*(C22**6)
        Q08=(C12**2)*(C22**7)
        Q09=C12*(C22**8)
        Q00=C22**9

        Q=np.matrix([[Q11, Q12, Q13, Q14, Q15, Q16, Q17, Q18, Q19, Q10], [Q21, Q22, Q23, Q24, Q25, Q26, Q27, Q28, Q29, Q20], [Q31, Q32, Q33, Q34, Q35, Q36, Q37, Q38, Q39, Q30], [Q41, Q42, Q43, Q44, Q45, Q46, Q47, Q48, Q49, Q40], [Q51, Q52, Q53, Q54, Q55, Q56, Q57, Q58, Q59, Q50], [Q61, Q62, Q63, Q64, Q65, Q66, Q67, Q68, Q69, Q60], [Q71, Q72, Q73, Q74, Q75, Q76, Q77, Q78, Q79, Q70], [Q81, Q82, Q83, Q84, Q85, Q86, Q87, Q88, Q89, Q80], [Q91, Q92, Q93, Q94, Q95, Q96, Q97, Q98, Q99, Q90], [Q01, Q02, Q03, Q04, Q05, Q06, Q07, Q08, Q09, Q00]], dtype=np.float)
        IQ=Q.I

        XQ1=A90*C11+B90*C12
        XQ2=A81*C11+B81*C12
        XQ3=A72*C11+B72*C12
        XQ4=A63*C11+B63*C12
        XQ5=A54*C11+B54*C12
        XQ6=A45*C11+B45*C12
        XQ7=A36*C11+B36*C12
        XQ8=A27*C11+B27*C12
        XQ9=A18*C11+B18*C12
        XQ0=A09*C11+B09*C12

        XQ=np.matrix([[XQ1], [XQ2], [XQ3], [XQ4], [XQ5], [XQ6], [XQ7], [XQ8], [XQ9], [XQ0]], dtype=np.float)
        MQ=np.dot(IQ, XQ)

        MQ1=MQ[0,0]
        MQ2=MQ[1,0]
        MQ3=MQ[2,0]
        MQ4=MQ[3,0]
        MQ5=MQ[4,0]
        MQ6=MQ[5,0]
        MQ7=MQ[6,0]
        MQ8=MQ[7,0]
        MQ9=MQ[8,0]
        MQ0=MQ[9,0]

        outhdr['PV1_49'] = np.float64(MQ1)
        outhdr['PV1_50'] = np.float64(MQ2)
        outhdr['PV1_51'] = np.float64(MQ3)
        outhdr['PV1_52'] = np.float64(MQ4)
        outhdr['PV1_53'] = np.float64(MQ5)
        outhdr['PV1_54'] = np.float64(MQ6)
        outhdr['PV1_55'] = np.float64(MQ7)
        outhdr['PV1_56'] = np.float64(MQ8)
        outhdr['PV1_57'] = np.float64(MQ9)
        outhdr['PV1_58'] = np.float64(MQ0)


        YQ1=A90*C21+B90*C22
        YQ2=A81*C21+B81*C22
        YQ3=A72*C21+B72*C22
        YQ4=A63*C21+B63*C22
        YQ5=A54*C21+B54*C22
        YQ6=A45*C21+B45*C22
        YQ7=A36*C21+B36*C22
        YQ8=A27*C21+B27*C22
        YQ9=A18*C21+B18*C22
        YQ0=A09*C21+B09*C22

        YQ=np.matrix([[YQ1], [YQ2], [YQ3], [YQ4], [YQ5], [YQ6], [YQ7], [YQ8], [YQ9], [YQ0]], dtype=np.float)
        NQ=np.dot(IQ, YQ)

        NQ1=NQ[9,0]
        NQ2=NQ[8,0]
        NQ3=NQ[7,0]
        NQ4=NQ[6,0]
        NQ5=NQ[5,0]
        NQ6=NQ[4,0]
        NQ7=NQ[3,0]
        NQ8=NQ[2,0]
        NQ9=NQ[1,0]
        NQ0=NQ[0,0]

        outhdr['PV2_49'] = np.float64(NQ1)
        outhdr['PV2_50'] = np.float64(NQ2)
        outhdr['PV2_51'] = np.float64(NQ3)
        outhdr['PV2_52'] = np.float64(NQ4)
        outhdr['PV2_53'] = np.float64(NQ5)
        outhdr['PV2_54'] = np.float64(NQ6)
        outhdr['PV2_55'] = np.float64(NQ7)
        outhdr['PV2_56'] = np.float64(NQ8)
        outhdr['PV2_57'] = np.float64(NQ9)
        outhdr['PV2_58'] = np.float64(NQ0)


    # r^9-terms
        outhdr['PV1_59'] = 0.
        outhdr['PV2_59'] = 0.

    return outhdr


def header_from_string(string):
    lines = [line for line in string.splitlines() 
             if len(line) == 80 and not line.startswith('#')]
    string = ''.join(lines)
    return fits.header.Header.fromstring(string)


def display_header(header):
    for kwd, value in header.items():
        print(f'{kwd:8s}= {value!s:>30s}')


def main():
    parser = argparse.ArgumentParser(
        description='Convert SIP to TPV non-linear coefficients')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='increase verbosity')
    parser.add_argument('-e', '--ext', default=0, type=int,
                        help='extension number')
    parser.add_argument('infile', nargs='?', default='-')
    parser.add_argument('outfile', nargs='?', default='-')
    args = parser.parse_args()
    if args.verbose == 0:
        loglevel = logging.WARNING
    elif args.verbose  == 1:
        loglevel = logging.INFO
    elif args.verbose >= 2:
        loglevel = logging.DEBUG
    logging.basicConfig(format='%(levelname)s:%(message)s', level=loglevel)

    logging.debug(f'Logging level set to {args.verbose}')
    logging.info('Reading SIP from %s', args.infile)
    logging.info('Writing TPV to %s', args.outfile)

    if args.infile == '-':
        string = sys.stdin.read()
    else:
        inhdr = fits.getheader(args.infile, args.ext)

    outhdr = sip2tpv(inhdr)

    primary = fits.PrimaryHDU()
    primary.header.update(outhdr)

    if args.outfile == '-':
        print(primary.header)
    else:
        primary.writeto(args.outfile)


if __name__ == '__main__':
    main()

