
'''
These functions are provided only to make fake data to test the workflow of the framework they are not validated and have no 
scientific validations to support the accuracy of the algorithms.

'''

def data_l2_density(conductivity, temp,pressure, lat, lon):
    '''
    Based on calculations done here: https://github.com/ooici/coi-services/blob/master/ion/processes/data/transforms/ctd/ctd_L2_density.py#L55
    '''
    from pygsw import vectors as gsw
    sp = gsw.sp_from_c(conductivity/42.914, temp, pressure)
    sa = gsw.sa_from_sp(sp, pressure, lon, lat)
    rho = gsw.rho(sa, temp, pressure)
    return rho

def data_l2_salinity(conductivity, temp, pressure):
    '''
    Based on calculations done here: https://github.com/ooici/coi-services/blob/master/ion/processes/data/transforms/ctd/ctd_L2_salinity.py#L58
    '''
    import pygsw.vectors as gsw
    sal_value = gsw.sp_from_c(conductivity/42.914, temp, pressure)
    return sal_value

