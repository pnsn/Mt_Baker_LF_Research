import datetime
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, SubplotSpec
import cartopy.crs as ccrs
from cartopy.io.img_tiles import GoogleTiles, OSM


WGS84 = ccrs.PlateCarree()
UTM10N = ccrs.UTM(zone=10, southern_hemisphere=False)


def rad2llur(rad=50000.):
    # Mount Baker Summit Coordinates
    lat= 48.7745 - 0.01
    lon=-121.8172 - 0.02

    # Convert reference location to northing & easting
    mEo, mNo = UTM10N.transform_point(lon,lat,WGS84)
    llE = mEo - rad
    llN = mNo - rad
    urE = mEo + rad
    urN = mNo + rad
    lowerleft = WGS84.transform_point(llE, llN, UTM10N)
    upperright = WGS84.transform_point(urE, urN, UTM10N)
    
    return [lowerleft[0], upperright[0], lowerleft[1], upperright[1]]

def mount_baker_basemap(
    fig=None, sps=None,
    radius_km=40, zoom=10,
    aws_terrain=True,
    open_street_map=True,
    aws_add_image_kwargs={'cmap':'Greys_r', 'alpha':0.5},
    osm_add_image_kwargs={}):

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
    extent = rad2llur(rad=radius_km*1e3)
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