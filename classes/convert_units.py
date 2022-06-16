
class conversion:
    """
        This class converts units from one unit system to another
    """
    def __init__(self,unit_from,unit_to):

        self.unitfrom = unit_from
        self.unitto = unit_to
        self.dtosi = {'si'     : self._si_to_si,
                      'cgs'    : self._cgs_to_si,
                      'micro'  : self._micro_to_si,
                      'nano'   : self._nano_to_si}

        self.dtounit = {'si'   : self._si_to_si,
                        'cgs'  : self._si_to_cgs,
                        'micro': self._si_to_micro,
                        'nano' : self._si_to_nano}

    @property
    def tosi(self):
        return self.dtosi[self.unitfrom]

    @property
    def tounit(self):
        return self.dtounit[self.unitto]

    def convert(self,property):
        return self.tosi[property] * self.tounit[property]

    @property
    def _si_to_si(self):
        return {'distance' :    1,      #m
                'mass':         1,      #kg
                'time':         1,      #s
                'energy':       1,      #kg m^2/s^2
                'velocity':     1,      #m/s
                'force':        1,      #kg m/s^2
                'torque':       1,      #kg m^2/s^2
                'dynviscosity': 1,      #kg/m
                'density':      1,      #kg/m^3
                'pressure':     1,      #kg/m*s
                'gamma' :       1,      #kg/m^2*s
                'viscosity':    1,      #kg/s
                'kn':           1,      #kg/m*s
                'tau':          1,      #kg/m*s
                'youngs':       1,      #kg/m*s
               }

    @property
    def _si_to_nano(self):
        return {'distance':     1e9,    #nm
                'mass':         1e21,   #ag
                'time':         1e9,    #ns
                'energy':       1e21,   #ag nm^2/ns^2
                'velocity':     1,      #nm/ns
                'force':        1e12,   #ag nm/ns^2
                'torque':       1e21,   #ag nm^2/ns^2
                'dynviscosity': 1e3,    #ag/nm*ns
                'density':      1e-6,   #ag/nm^3
                'pressure':     1e-6,   #ag/nm*ns^2
                'gamma' :       1e-6,   #ag/nm^2*ns
                'viscosity':    1e12,   #ag/ns
                'kn':           1e-6,   #ag/nm*ns^2
                'tau':          1e-6,   #ag/nm*ns^2
                'youngs':       1e-6,   #ag/nm*ns^2
               }

    @property
    def _si_to_cgs(self):
        return {'distance' :    1e2,    #cm
                'mass':         1e3,    #g
                'time':         1,      #s
                'energy':       1e7,    #ergs
                'velocity':     1e2,    #cm/s
                'force':        1e5,    #dynes
                'torque':       1e7,    #dyne*centimeters
                'dynviscosity': 1e1,    #Poise
                'density':      1e-3,   #gram/cm^3
                'pressure':     1e1,    #dyne/cm^2
                'gamma' :       1e-1,   #g/cm^2*s
                'viscosity':    1e3,    #g/s
                'kn':           1e1,    #dyne/cm^2
                'tau':          1e1,    #dyne/cm^2
                'youngs':       1e1,    #dyne/cm^2
               }

    @property
    def _si_to_micro(self):
        # u=micro
        return {'distance' :    1e6,    #um
                'mass':         1e15,   #pg
                'time':         1e6,    #us
                'energy':       1e15,   #pg*um^2/us^2
                'velocity':     1,      #um/us
                'force':        1e9,    #pg*um/us^2
                'torque':       1e15,   #pg*um^2/us^2
                'dynviscosity': 1e3,    #pg/um*us
                'density':      1e-3,   #pg/um^2
                'pressure':     1e-3,   #pg/um*us^2
                'gamma' :       1e-2,   #pg/um^2*us
                'viscosity':    1e9,    #pg/us
                'kn':           1e-3,   #pg/um*us^2
                'tau':          1e-3,   #pg/um*us^2
                'youngs':       1e-3,   #pg/um*us^2
               }

    @property
    def _nano_to_si(self):
        return {'distance':     1e-9,   #nm
                'mass':         1e-21,  #ag
                'time':         1e-9,   #ns
                'energy':       1e-21,  #ag nm^2/ns^2
                'velocity':     1,      #nm/ns
                'force':        1e-12,  #ag nm/ns^2
                'torque':       1e-21,  #ag nm^2/ns^2
                'dynviscosity': 1e-3,   #ag/nm*ns
                'density':      1e6,    #ag/nm^3
                'pressure':     1e6,    #ag/nm*ns^2
                'gamma' :       1e6,    #ag/nm^2*ns
                'viscosity':    1e-12,  #ag/ns
                'kn':           1e6,    #ag/nm*ns^2
                'tau':          1e6,    #ag/nm*ns^2
                'youngs':       1e6,    #ag/nm*ns^2
               }

    @property
    def _cgs_to_si(self):
        return {'distance' :    1e-2,   #cm
                'mass':         1e-3,   #g
                'time':         1,      #s
                'energy':       1e-7,   #ergs
                'velocity':     1e-2,   #cm/s
                'force':        1e-5,   #dynes
                'torque':       1e-7,   #dyne*centimeters
                'dynviscosity': 1e-1,   #Poise
                'density':      1e3,    #gram/cm^3
                'pressure':     1e-1,   #dyne/cm^2
                'gamma' :       1e1,    #g/cm^2*s
                'viscosity':    1e-3,   #g/s
                'kn':           1e-1,   #dyne/cm^2
                'tau':          1e-1,   #dyne/cm^2
                'youngs':       1e-1,   #dyne/cm^2
               }

    @property
    def _micro_to_si(self):
        # u=micro
        return {'distance' :    1e-6,   #um
                'mass':         1e-15,  #pg
                'time':         1e-6,   #us
                'energy':       1e-15,  #pg*um^2/us^2
                'velocity':     1,      #um/us
                'force':        1e-9,   #pg*um/us^2
                'torque':       1e-15,  #pg*um^2/us^2
                'dynviscosity': 1e-3,   #pg/um*us
                'density':      1e3,    #pg/um^2
                'pressure':     1e3,    #pg/um*us^2
                'gamma' :       1e2,    #pg/um^2*us
                'viscosity':    1e-9,   #pg/us
                'kn':           1e3,    #pg/um*us^2
                'tau':          1e3,    #pg/um*us^2
                'youngs':       1e3,    #pg/um*us^2
               }
