# -*- coding: utf-8 -*-
"""
Ethernet_connection_by_MAC.py 31.01.2022 Rev0

This example script tries to establish a connection to a NV200D/NET when only 
the MAC address is known. The IP address is determined by sending a UDP 
broadcast and receiving answers form Lantronix modules. 


"""

from socket import socket, AF_INET,SOCK_DGRAM,IPPROTO_UDP,SOL_SOCKET,SO_BROADCAST # required for UDP broadcast
import telnetlib  # required for Telnet communication

# define communication parameters
#MAC = '00:80:A3:EC:B1:7A' # desired MAC address of NV200D/NET
MAC = '00:80:A3:E5:C4:E6'
port=23                   # Telnet Port (default: 23)
timeout=0.4               # Timeout for reading device answers


data = None     # will be assigned with UDP broadcast response
IP = None       # will be assigned with the IP address
NV200 = None    # will be assigned with the Telnet handle

# Send a UDP broadcast to network to get a response from all NV200D/NET devices
print('Sending UDP broadcast')
s=socket(AF_INET, SOCK_DGRAM,IPPROTO_UDP)
s.setsockopt(SOL_SOCKET, SO_BROADCAST, 1)
s.settimeout(timeout)
s.sendto(bytearray([0x00,0x00,0x00,0xf6]),('255.255.255.255',30718)) # using Lantronix Discovery Protocol 

data=[]
try:
    while(True): 
        recvbytes = s.recvfrom(65565)   # receive message (msg, src address(ip, port))
        data.append(recvbytes)
except:
    if not data:  
        print('Timeout occured, no response to broadcast')

# Print a list of all responding devices IP and MAC addresses
if data:
    devices=[]
    for da in data:       
      helper=str(da[0].hex())[48:].upper()  # mac address (the last 12 byte)
      helper=':'.join(helper[i:i+2] for i in range(0,12,2))  
      devices.append({'MAC':helper,'IP': str(da[1][0])})
    # Print IP and MAC addresses of all devices found 
    print('%i device(s) responded:'%len(devices))
    for dev in devices:
        print('\t%s\t%s'%(dev['IP'],dev['MAC']))
        if dev['MAC']==MAC: # assign the IP address for the desired MAC address
            IP=dev['IP']

# Establish Telnet connection by IP address
if IP:
    NV200=telnetlib.Telnet(IP,port)
    print('Established Telnet connection to IP: %s'%IP)
else:
    print('Desired device MAC address %s not in list'%MAC)

# Query, read and print device prompt, then close connection 
if NV200:
    NV200.write(b'\r')   
    print(NV200.read_until(b'\n',timeout).decode(),end='')
    NV200.close()

