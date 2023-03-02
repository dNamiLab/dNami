
import numpy as np
import rasterio as rio

fbathy = 'bathy.tif'
try:
    with rio.open(fbathy) as f:
        bathy = f.read(1)
except Exception as e:
    print(e)
    exit()
bathy = np.flipud(bathy)

with open('bathy.bin', 'wb') as f:
    np.save(f, bathy)
