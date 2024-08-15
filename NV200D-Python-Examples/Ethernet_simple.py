# -*- coding: utf-8 -*-
"""
Ethernet_simple.py 23.12.2021 Rev0

This example script shows the basic way to communicate with the NV200D/NET via 
Telnet with known IP address. First the connection is established and it is checked if the 
controller reacts. Then a closed loop step to a position of 10 µm (or 10mrad 
for tilting actuators) is executed and after one second of waiting time the
reached position of the actuator is queried and displayed.

Expected output:
NV200/D_NET>
meas,9.942

Note that the value may vary due to sensor noise and controller performance

"""


import telnetlib # required for Telnet communication
import time      # required to wait for one second


IP='192.168.188.71'     # NV200D/NET IP address, has to be adjusted to your address
port=23                 # Telnet Port (default: 23)
timeout=0.25            # Timeout for reading device answers

NV200 = telnetlib.Telnet(IP, port)  # open Telnet connection to NV200D/NET

NV200.write(b'\r')   # Query device prompt response 
print(NV200.read_until(b'\n',timeout).decode(),end='') # Read and print device prompt response

NV200.write(b'modsrc,0\r')   # Set setpoint source to Telnet/USB input
NV200.write(b'cl,1\r')       # Set controller to closed loop mode
NV200.write(b'set,10\r')     # Set setpoint to 10 µm (or mrad for tilting actuators)
time.sleep(1)                               # Wait for actuator to change position
NV200.write(b'meas\r')       # Query position measurement
print(NV200.read_until(b'\n',timeout).decode(),end='') # Read and print device answer to position measurement query

NV200.close()   # close Telnet connection to NV200D/NET