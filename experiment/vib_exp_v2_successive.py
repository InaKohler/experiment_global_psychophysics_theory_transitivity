# -*- coding: utf-8 -*-
"""
vib_exp_v2_successive.py

Connects to the tactor controller (Engineering Acoustics, C-2 tactor) via
TactorInterface.dll on COM6.

The module is used only for its side effect: importing it (in
run_experiment.py, global_setup.py and trialFunction.py) runs the connection
code below once. The other modules load the same DLL and send the vibration
pulses themselves.

Created:  2025/02/26, Marcus Dietlein
Last mod: 2026/09/30, Ina Kohler (reduced to the tactor connection)
"""

import time
from ctypes import *

from psychopy import prefs
prefs.hardware['audioLib'] = ['PTB']        # use psychtoolbox (PTB) for python3
prefs.general['shutdownKey'] = 'q'

# Load library
lib = CDLL('TactorInterface.dll')

# Initialize TI
lib.InitializeTI()
print('Initialized TI')

# Add a delay after initialization
time.sleep(2)  # Give the system time to stabilize

# Connect to the controller on COM6 in USB mode
print("Attempting connection to COM6...")
try:
    lib.Discover(c_int(1))
    time.sleep(1)  # Wait after discovery

    port = "COM6"
    encoded_port = port.encode('utf-8')
    resp = lib.Connect(c_char_p(encoded_port), c_int(1), None)
    print(f"Connection response: {resp}")

    if resp >= 0:
        print("Successfully connected!")
    else:
        print("Connection failed")
        lib.ShutdownTI()
        raise RuntimeError(f"Could not connect to device. Error code: {resp}")
except Exception as e:
    lib.ShutdownTI()
    raise
