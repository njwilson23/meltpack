import karta
import numpy as np

def overlap_bbox(bbox1, bbox2):
    bbox = [max(bbox1[0], bbox2[0]), max(bbox1[1], bbox2[1]),
            min(bbox1[2], bbox2[2]), min(bbox1[3], bbox2[3])]
    if (bbox[0] > bbox[2]) or (bbox[1] > bbox[3]):
        raise ValueError("input bounding boxes not overlapping")
    return bbox

def apply_shared_pixels(func, scene1, scene2, nodata=None):
    """ Apply a function *f(a,b) -> c* to pixels shared between two scenes. """
    extent = scene1.get_data_extent()
    bbox1 = (extent[0], extent[2], extent[1], extent[3])
    extent = scene2.get_data_extent()
    bbox2 = (extent[0], extent[2], extent[1], extent[3])

    if nodata is None:
        nodata = scene1.nodata

    if np.isnan(nodata):
        isnodata = np.isnan
    else:
        def isnodata(a):
            return a == nodata

    assert(isnodata(scene2.nodata))         # require scenes to use same nodata

    bbox = overlap_bbox(bbox1, bbox2)
    dx, dy = scene1.transform[2:4]
    bbox = (bbox[0]+dx/2, bbox[1]+dy/2, bbox[2]-dx/2, bbox[3]-dy/2)

    # Calculate the bounds of data pixels for each image
    idx_bounds1 = (scene1.get_indices(bbox[0], bbox[1]),
                   scene1.get_indices(bbox[2], bbox[3]))
    idx_bounds2 = (scene2.get_indices(bbox[0], bbox[1]),
                   scene2.get_indices(bbox[2], bbox[3]))

    # Choose an appropriate number of rows to read as a time
    nx = idx_bounds1[1][1] - idx_bounds1[0][1]
    ny = idx_bounds1[1][0] - idx_bounds1[0][0]

    if nx != (idx_bounds2[1][1] - idx_bounds2[0][1]):
        print("WARNING: incompatible swath sizes in count_shared_pixels")

    bandwidth = 1000000 // nx      # assuming 4-byte numbers and 400 mb memory use (conservative)
    output = np.nan*np.zeros([ny, nx], dtype=scene1.values.dtype)

    i = 0
    while i != ny:
        inext = min(ny, i+bandwidth)

        swath1 = scene1.values[idx_bounds1[0][0]+i:idx_bounds1[0][0]+inext,
                               idx_bounds1[0][1]:idx_bounds1[1][1]]

        swath2 = scene2.values[idx_bounds2[0][0]+i:idx_bounds2[0][0]+inext,
                               idx_bounds2[0][1]:idx_bounds2[1][1]]

        output[i:inext,:] = func(swath1, swath2)
        i = inext

    return output[~np.isnan(output)]

def difference_shared_pixels(scene1, scene2, nodata=None):
    """ Convenience function to return differences between pixels shared
    between two grids. """
    return apply_shared_pixels(lambda a, b: a-b, scene1, scene2, nodata=nodata)

def count_shared_pixels(scene1, scene2, bbox1=None, bbox2=None, nodata=None):
    """ Count the valid data pixels shared between two scenes.

    Arguments:
    ----------
    scene1: RegularGrid
    scene2: RegularGrid
    bbox1: tuple, optional
    bbox2: tuple, optional
    nodata: nodata value, optional

    Returns:
    --------
    int
    """
    if bbox1 is None:
        extent = scene1.get_data_extent()
        bbox1 = (extent[0], extent[2], extent[1], extent[3])
    if bbox2 is None:
        extent = scene2.get_data_extent()
        bbox2 = (extent[0], extent[2], extent[1], extent[3])
    if nodata is None:
        nodata = scene1.nodata

    if np.isnan(nodata):
        isnodata = np.isnan
    else:
        def isnodata(a):
            return a == nodata

    assert(isnodata(scene2.nodata))         # require scenes to use same nodata

    try:
        bbox = overlap_bbox(bbox1, bbox2)
        dx, dy = scene1.transform[2:4]
        bbox = (bbox[0]+dx/2, bbox[1]+dy/2, bbox[2]-dx/2, bbox[3]-dy/2)
    except ValueError:
        return 0

    # Calculate the bounds of data pixels for each image
    idx_bounds1 = (scene1.get_indices(bbox[0], bbox[1]),
                   scene1.get_indices(bbox[2], bbox[3]))
    idx_bounds2 = (scene2.get_indices(bbox[0], bbox[1]),
                   scene2.get_indices(bbox[2], bbox[3]))

    # Choose an appropriate number of rows to read as a time
    nx = idx_bounds1[1][1] - idx_bounds1[0][1]
    ny = idx_bounds1[1][0] - idx_bounds1[0][0]

    if nx != (idx_bounds2[1][1] - idx_bounds2[0][1]):
        print("WARNING: incompatible swath sizes in count_shared_pixels")

    bandwidth = 1000000 // nx      # assuming 4-byte numbers and 400 mb memory use (conservative)
    pixel_count = 0
    i = 0
    while i != ny:
        inext = min(ny, i+bandwidth)

        swath1 = scene1.values[idx_bounds1[0][0]+i:idx_bounds1[0][0]+inext,
                               idx_bounds1[0][1]:idx_bounds1[1][1]]

        swath2 = scene2.values[idx_bounds2[0][0]+i:idx_bounds2[0][0]+inext,
                               idx_bounds2[0][1]:idx_bounds2[1][1]]
        i = inext

        pixel_count += np.sum((~isnodata(swath1)) & (~isnodata(swath2)))

    return pixel_count

