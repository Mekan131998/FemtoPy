# -*- coding: utf-8 -*-
"""
USB_simple.py 23.12.2021 Rev0

This example script shows the basic way to communicate with the NV200D/NET via 
USB. First the connection is established and it is checked if the 
controller reacts. Then a closed loop step to a position of 10 µm (or 10mrad 
for tilting actuators) is executed and after one second of waiting time the
reached position of the actuator is queried and displayed.

Expected output:
NV200/D_NET>
meas,9.942

Note that the value may vary due to sensor noise and controller performance

"""

import serial    # required for serial communication, requires "pip install pyserial"
import time      # required to wait for one second


COM='COM10'         # NV200D/NET COM port address, has to be adjusted to your address
baudrate=115200     # Baudrate (default: 115200)
timeout=0.25        # Timeout for reading device answers

NV200 = serial.Serial(COM, baudrate, timeout=timeout, xonxoff=True) # open serial connection to NV200D/NET


NV200.write(b'\r')   # Query device prompt response 
print(NV200.readline().decode(),end='') # Read and print device prompt response

NV200.write(b'modsrc,0\r')   # Set setpoint source to Telnet/USB input
NV200.write(b'cl,1\r')       # Set controller to closed loop mode
NV200.write(b'set,10\r')     # Set setpoint to 10 µm (or mrad for tilting actuators)
time.sleep(1)                               # Wait for actuator to change position
NV200.write(b'meas\r')       # Query position measurement
print(NV200.readline().decode(),end='')     # Read and print device answer to position measurement query

NV200.close()   # close Telnet connection to NV200D/NET