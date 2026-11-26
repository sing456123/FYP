import socket
s = socket.create_connection(("ntrip.geodetic.gov.hk", 2101), timeout=10)
s.sendall(b"GET /HKSC_32 HTTP/1.0\r\nHost: ntrip.geodetic.gov.hk\r\n"
          b"Ntrip-Version: Ntrip/2.0\r\nUser-Agent: NTRIP py\r\n\r\n")
hdr = b""
while b"\r\n\r\n" not in hdr:
    hdr += s.recv(1024)
print(hdr)                       
data = s.recv(4096)

print(data)