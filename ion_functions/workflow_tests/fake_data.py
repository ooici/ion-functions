
'''
These functions are provided only to make fake data to test the workflow of the framework they are not validated and have no 
scientific validations to support the accuracy of the algorithms.

'''

def data_l2_density(conductivity, temp,pressure, lat, lon):
    '''
    '''
    from pygsw import vectors as gsw
    # 
    sp = gsw.sp_from_c(conductivity, temp, pressure)
    sa = gsw.sa_from_sp(sp, pressure, lon, lat)
    rho = gsw.rho(sa, temp, pressure)
    return rho

def data_l2_salinity(conductivity, temp, pressure):
    '''
    '''
    import pygsw.vectors as gsw
    sal_value = gsw.sp_from_c(conductivity, temp, pressure)
    return sal_value

