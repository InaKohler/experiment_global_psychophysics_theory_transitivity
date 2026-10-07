"""
conversions.py

Conversions between the intensity units of the visual stimuli using the
monitor's lookup table:
  visualscale : device units used in session files and during adjustment
  rgb         : grey value sent to the monitor (0-255)
  db          : luminance in dB Lambert (measured)

The lookup table lookuptable/lookuptab_mavo.txt has one row per RGB value
with the columns rgb, db and visualscale (separator ";").

Last mod: 30/09/2026, Ina Kohler
"""
from numpy import argmin
from pandas import read_table

with open("lookuptable/lookuptab_mavo.txt", encoding = "utf8") as lookupfile:
    lookuptab = read_table(lookupfile, sep = ";")

def db_to_rgb(intensity_db):
    """Return the RGB value whose luminance is closest to intensity_db. The
    row index of the lookup table is identical with the RGB value.
    """
    return argmin(abs(lookuptab["db"] - intensity_db))

def visualscale_to_rgb(intensity_visualscale):
    """Return the RGB value for a visualscale intensity.
    """
    return lookuptab.loc[(lookuptab['visualscale'] == intensity_visualscale),
                         "rgb"].values[0]

def visualscale_to_db(intensity_visualscale):
    """Return the luminance in dB Lambert for a visualscale intensity.
    """
    return lookuptab.loc[(lookuptab['visualscale'] == intensity_visualscale),
                         "db"].values[0]
