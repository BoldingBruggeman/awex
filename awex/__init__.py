import enum

from ._awex import *


class HumidityMeasure(enum.IntEnum):
    """Measure used to specify air humidity"""

    RELATIVE_HUMIDITY = 1  #: relative humidity in %
    WET_BULB_TEMPERATURE = 2  #: wet bulb temperature in in degrees Celsius
    DEW_POINT_TEMPERATURE = 3  #: dewpoint temperature in in degrees Celsius
    SPECIFIC_HUMIDITY = 4  #: specific humidity in kg kg-1


class LongwaveMethod(enum.IntEnum):
    """Method used to calculate longwave radiation"""

    NET = 1  #: Surface net longwave radiation - ERA5 - str
    DOWNWARDS = 2  #: Surface downwards longwave radiation downwards - ERA5 - strd
    CLARK = 3  #: Clark et al. (1974)
    HASTENRATH_LAMB = 4  #: Hastenrath and Lamb (1978)
    BIGNAMI = 5  #: Bignami et al. (1995)
    BERLIAND_BERLIAND = 6  #: Berliand & Berliand (1952)
    JOSEY1 = 7  #: Josey et.al. 2003 - (J1,9)
    JOSEY2 = 8  #: Josey et.al. 2003 - (J2,14)


class ShortwaveMethod(enum.IntEnum):
    """Method used to calculate shortwave radiation"""

    NET = 1  #: Surface net shortwave radiation - ERA5 - ssr
    DOWNWARDS = 2  #: Surface shortwave radiation downwards - ERA5 - ssrd
    ROSATI_MIYAKODA = 3  #: Rosati & Miyakodi (1988)


class AlbedoMethod(enum.IntEnum):
    """Method used to calculate albedo"""

    PAYNE = 1  #: Payne (1972)
    COGLEY = 2  #: Cogley (1979)
