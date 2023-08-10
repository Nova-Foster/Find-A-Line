import numpy as np
import codecs

#raw_imd = open("800V - 250Hz - 500usmm - 15um vertical.imd","r",encoding="hex_codec")
raw = codecs.open("800V - 250Hz - 500usmm - 15um vertical.imd","r",encoding="utf16")
print(raw)

test = raw.readline()
print(test)