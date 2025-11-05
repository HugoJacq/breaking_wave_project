import xarray as xr

# read netcdf, set up for xgcm
def read_data(filename: string):
    """
    This function reads the output of Basilisk (using bderembl/libs/netcdf_pas.h)
        and return a dataset with xgcm coordinates to allow for easy differentiation

    INPUT:
        filename: a string where the .nc is
    OUPUT:
        dataset: a xr.Dataset
    """
    ds_out = xr.open_mfdataset(filename)
    
    # building a new dataset with xgcm grid


