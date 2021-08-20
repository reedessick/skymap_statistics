description = "a module that houses a set of colors"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

colors = 'b g r c m k y'.split()
N = len(colors)

ifoColors = {'H':'r',
             'L':'b',
             'V':'g',
             'G':'k',
             'K':'m',
            }

#-------------------------------------------------

def getColor():
    '''
    iterates through colors in a set order
    '''
    i = 0
    while 1:
        yield colors[i]
        i = (i+1)%N

def getIFOColor( ifo ):
    ''' 
    returns a specific color for each known IFO.
    if IFO is not know, the default is 'k'
    '''
    if ifo in ifoColors:
        return ifoColors[ifo]
    else:
        return 'k'
