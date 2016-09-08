description = "a module that houses classes that write html structures for fits2html.py and friends"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

class snglFITS(dict):
    '''
    a class that houses data and renders html, json for info about a single FITS file (ie: no comparisons)
    '''

    def __init__(self, fitsname, *args):
        self.fitsname = fitsname
        super(snglFITS, self).__init__(*args)

class multFITS(dict):
    '''
    a class that houses data and renders html, json for comparison info about multiple FITS files
    '''

    def __init__(self, *args):
        self.fitsnames = []
        super(multFITS, self).__init__(*args)
