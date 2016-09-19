description = "a module that houses a set of colors"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

colors = 'b g r c m k y'.split()
N = len(colors)

#-------------------------------------------------

def getColor():
    i = 0
    while 1:
        yield colors[i]
        i = (i+1)%N
