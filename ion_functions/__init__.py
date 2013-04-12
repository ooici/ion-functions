from version import version as __version__

from utils import isempty, islogical, ismatrix, isnumeric, isreal, isscalar, isvector
from qc.qc_functions import dataqc_condcompress, dataqc_globalrangetest, dataqc_globalrangetest_minmax, dataqc_gradienttest, dataqc_localrangetest, dataqc_polytrendtest, dataqc_propogateflags, dataqc_solarelevation, dataqc_spiketest, dataqc_stuckvaluetest
from data.ctd_functions import ctd_density, ctd_pracsal
from data.adcp_functions import adcp_beam2ins, adcp_beam_eastward, adcp_beam_error, adcp_beam_northward, adcp_beam_vertical, adcp_ins2earth, adcp_magvar

