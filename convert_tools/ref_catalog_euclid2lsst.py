#! /usr/bin/env python3
'''
Created on Mar 11, 2021

@author: colley
'''
import os
import sys
import astropy.io.fits as fits
from astropy import table
import numpy as np
from astropy.table import Column


def convert_jy_mag(t_tab):
    '''Convert Jy to magnitude
    Le 11/03/2021 a 12:17, Santiago Serrano Elorduy
    You can easily recover the magnitude from the FNU fluxes with 
        m = -2.5*log10(TU_FNU_XXX / 3631)
    :param t_tab: tab
    '''    
    return -2.5 * np.log10(t_tab / 3631.0)
    
    
def ref_cat_conversion(input_cat, output_cat_file):
    '''read, convert and write reference cataloog for stack LSST
    
    :param input_cat:
    :param output_cat_file:
    '''
    with fits.open(input_cat) as hdul:
        cat_table = hdul[1].data
    nb_star = len(cat_table)
    print("nb_star: ", nb_star)
    o_table = table.Table()
    o_table.add_column(table.Column(data=range(nb_star), name="uniqueId"))
    # cop
#     print(cat_table["RA"])
#     print(type(cat_table["RA"]))
    o_table.add_column(Column(cat_table["RA"]), name="RA")
    o_table.add_column(Column(cat_table["DEC"]), name="DEC")
    o_table.add_column(Column(cat_table["TU_FNU_U_LSST"]), name="u")
    o_table.add_column(Column(cat_table["TU_FNU_G_LSST"]), name="g")
    o_table.add_column(Column(cat_table["TU_FNU_R_LSST"]), name="r")
    o_table.add_column(Column(cat_table["TU_FNU_I_LSST"]), name="i")
    o_table.add_column(Column(cat_table["TU_FNU_Z_LSST"]), name="z")
    o_table.add_column(Column(cat_table["TU_FNU_Y_LSST"]), name="y")
    # print(o_table)
    del cat_table
    # convert Jy to mag
    o_table["u"] = convert_jy_mag(o_table["u"])    
    o_table["g"] = convert_jy_mag(o_table["g"])    
    o_table["r"] = convert_jy_mag(o_table["r"])    
    o_table["i"] = convert_jy_mag(o_table["i"])    
    o_table["z"] = convert_jy_mag(o_table["z"])    
    o_table["y"] = convert_jy_mag(o_table["y"])    
    # o_table.add_column(table.Column(data=np.ones(nb_star, dtype=int), name="isresolved")) 
    o_table.add_column(table.Column(data=np.zeros(nb_star, dtype=int), name="isvariable"))
    # print(o_table)
    # o_table.write(output_cat_file, format="ascii.commented_header", delimiter=',')
    return o_table, nb_star
    
    
def ref_cat_conversion_2018(input_cat, output_cat_file):
    o_table, nb_star = ref_cat_conversion(input_cat, output_cat_file)
    o_table.add_column(table.Column(data=np.ones(nb_star, dtype=int), name="isresolved"))
    o_table.write(output_cat_file, format="ascii.commented_header", delimiter=',')


def ref_cat_conversion_2021(input_cat, output_cat_file):
    o_table, nb_star = ref_cat_conversion(input_cat, output_cat_file)
    o_table.rename_column("y", "Y")
    o_table.add_column(table.Column(data=np.zeros(nb_star, dtype=int), name="resolved"))
    o_table.write(output_cat_file, format="ascii.commented_header", delimiter=',')

    
if __name__ == "__main__":
    input_cat_file = sys.argv[1]
    output_cat_file = sys.argv[2]
    os.system('rm -rf ' + output_cat_file)
    ref_cat_conversion(input_cat_file, output_cat_file)
