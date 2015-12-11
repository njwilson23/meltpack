import karta
import numpy as np

def count_shared_pixels(scene1, scene2, bbox1=None, bbox2=None, nodata=np.nan):
    if bbox1 is None:
        extent = scene1.get_data_extent()
        bbox1 = (extent[0], extent[2], extent[1], extent[3])
    if bbox2 is None:
        extent = scene2.get_data_extent()
        bbox2 = (extent[0], extent[2], extent[1], extent[3])

    nodata = scene1.nodata
    if np.isnan(nodata):
        isnodata = np.isnan
    else:
        def isnodata(a):
            return a == nodata

    assert(isnodata(scene2.nodata))         # require scenes to use same nodata

    # Calculate the bbox of overlap
    bbox = [max(bbox1[0], bbox2[0]), max(bbox1[1], bbox2[1]),
            min(bbox1[2], bbox2[2]), min(bbox1[3], bbox2[3])]

    if bbox[0] > bbox[2]:
        return 0
    if bbox[1] > bbox[3]:
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

