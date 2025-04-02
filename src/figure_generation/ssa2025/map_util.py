import datetime
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, SubplotSpec
import matplotlib.colors as colors
import cartopy.crs as ccrs
from cartopy.io.img_tiles import GoogleTiles, OSM
import cartopy.feature as cfeature

WGS84 = ccrs.PlateCarree()
UTM10N = ccrs.UTM(zone=10, southern_hemisphere=False)

BAKER_LAT = 48.7745
BAKER_LON = -121.8172

def rad2llur(rad=50000., latnudge=0, lonnudge=0):
    # Mount Baker Summit Coordinates
    lat = BAKER_LAT + latnudge
    lon = BAKER_LON + lonnudge

    # Convert reference location to northing & easting
    mEo, mNo = UTM10N.transform_point(lon,lat,WGS84)
    llE = mEo - rad
    llN = mNo - rad
    urE = mEo + rad
    urN = mNo + rad
    lowerleft = WGS84.transform_point(llE, llN, UTM10N)
    upperright = WGS84.transform_point(urE, urN, UTM10N)
    
    return [lowerleft[0], upperright[0], lowerleft[1], upperright[1]]


def radiusllsets(rad=50000., npts=101):
    lat = BAKER_LAT
    lon = BAKER_LON
    mEo, mNo = UTM10N.transform_point(lon, lat, WGS84)
    dtheta = 2.*np.pi/(npts - 1)
    theta = np.linspace(0, 2.*np.pi + dtheta, npts)
    mE = np.cos(theta)*rad + mEo
    mN = np.sin(theta)*rad + mNo
    return (mE, mN)


def mount_baker_basemap(
    fig=None, sps=None,
    radius_km=40, zoom=10,
    aws_terrain=True,
    open_street_map=True,
    aws_add_image_kwargs={'cmap':'Greys_r', 'alpha':0.5},
    osm_add_image_kwargs={},
    latnudge=-0.01,lonnudge=-0.02):

    if not aws_terrain and not open_street_map:
        raise ValueError

    # Establish Tile Query Services
    # Amazon Web Services elevation tiles
    attrib = []
    if aws_terrain:
        aws = GoogleTiles(
            desired_tile_form='L',
            url=('https://s3.amazonaws.com/elevation-tiles-prod/normal/{z}/{x}/{y}.png'),
        )
        
        _attrib = f" Terrain Tiles was accesed on {datetime.datetime.now().strftime('%Y-%m-%d')}"
        _attrib += '\n    from https://registry.opendata.aws/terrain-tiles'
        attrib.append(_attrib)
        imagery = aws

    if open_street_map:
        osm = OSM(cache=True)
        _attrib = f" Open Street Map accessed on {datetime.datetime.now().strftime('%Y-%m-%d')} using CartoPy"
        attrib.append(_attrib)
        if not aws_terrain:
            imagery = osm

    if fig is None:
        fig = plt.figure()

    if sps is None:
        ax = fig.add_subplot(111, projection=imagery.crs)
    elif isinstance(sps, SubplotSpec):
        ax = fig.add_subplot(sps, projection=imagery.crs)
    else:
        raise TypeError


    # Set map extent & get imagery
    extent = rad2llur(rad=radius_km*1e3, latnudge=latnudge, lonnudge=lonnudge)
    ax.set_extent(extent, ccrs.Geodetic())

    if open_street_map:
        ax.add_image(osm, zoom, **osm_add_image_kwargs)
        # print('added OSM')

    if aws_terrain:
        ax.add_image(aws, zoom, **aws_add_image_kwargs)
        # print('added AWS')
    


    if len(attrib) > 0:
        return (ax, attrib)
    else:
        return (ax, None)
    

def add_inset_map(fig, extent=[0.725, 0.7, 0.15, 0.15], pad_deg=10, projection=ccrs.PlateCarree()):
    axs = fig.add_axes(extent, projection=projection)
    axs.set_extent([BAKER_LON - pad_deg, BAKER_LON + pad_deg, 
                    BAKER_LAT - pad_deg, BAKER_LAT + pad_deg], 
                    crs=projection)
    axs.add_feature(cfeature.LAND)
    axs.add_feature(cfeature.OCEAN)
    axs.add_feature(cfeature.COASTLINE)
    axs.add_feature(cfeature.BORDERS, linestyle='-')
    axs.add_feature(cfeature.STATES, linestyle=':')
    axs.add_feature(cfeature.LAKES, alpha=0.5)
    return axs


# def add_grids(geoaxis, del_lon=0.25, del_la)

def add_rings(
        geoaxis,
        rads_km=[10, 30, 50, 70, 90],
        rads_colors=['black','firebrick','darkgreen','darkblue','purple'],
        npts=101,
        label_pt=28,
        annotations=[None, None, None, None, None],
        **options):
    for _e, _r in enumerate(rads_km):
        # Plot circles
        mE, mN = radiusllsets(rad=_r*1e3, npts=npts)
        geoaxis.plot(mE, mN, ':', color=rads_colors[_e], alpha=0.667, transform=UTM10N)
        # Label Radii
        text = f'{_r:d} km'
        if 'va' not in options.keys():
            options.update({'va':'bottom'})
        if 'ha' not in options.keys():
            options.update({'ha': 'left'})
        if not annotations[_e] is None:
            text += f'\n{annotations[_e]}'
        geoaxis.text(mE[label_pt], mN[label_pt]+1e3,text, transform=UTM10N, **options)



def mark_mount_baker(geoaxis, labeled=True, xoffset=0.05, yoffset=0, zorder=1000):
    kw = {'marker': '^', 'ms': 6, 'mec': 'k', 'mfc': 'w', 'transform': ccrs.PlateCarree(), 'zorder': zorder}
    geoaxis.plot(BAKER_LON, BAKER_LAT, **kw)
    if labeled:
        geoaxis.text(BAKER_LON + xoffset, BAKER_LAT + yoffset, 'Mount\nBaker', va='top', transform=ccrs.PlateCarree())

def plot_baker(geoaxis, zorder=2):
    handle = geoaxis.scatter(
        BAKER_LON, BAKER_LAT, marker='^',
        s=64, c='orange',
        edgecolors='k', zorder=zorder, 
        label='Mount Baker',
        transform=ccrs.PlateCarree())
    return handle

def pnsn_pallet():
    pp = {'evergreen': '#107a10',
          'black': 'k',
          'mid gray': '#ABABAB',
          'light gray': '#F5F5F5',
          'white': 'w',
          'forest green': '#094309',
          'mint': '#CAE3CA',
          'light blue': '#BCE6F8',
          'lime': '#53D623',
          'navy': '#00425D'}
    return pp

def make_pnsn_cmap(
        pallet_names=['lime','evergreen','forestgreen'],
        discretization=None,
        cmap_handle='pnsn_greens'):
    color_list = [pnsn_pallet()[_k] for _k in pallet_names]
    if discretization is None:
        mycmap = colors.LinearSegmentedColormap.from_list(
            cmap_handle, color_list)
        return mycmap
    elif len(color_list) == len(discretization) - 1:
        mycmap = colors.ListedColormap(color_list)
        norm = colors.BoundaryNorm(discretization, mycmap.N)
        return mycmap, norm

