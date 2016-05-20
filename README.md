# QSEIS (packaged as Fomosto backend)

Code to calculate synthetic seismograms based on a layered viscoelastic
half-space earth model.

QSEIS2D 2014 has been written by Rongjiang Wang.

Packaging has been done by Sebastian Heimann.

## References

Wang, R., (1999), A simple orthonormalization method for stable and efficient
computation of Green's functions, Bulletin of the Seismological Society of
America, 89(3), 733-741.

## Compile and install

```
autoreconf -i   # only if 'configure' script is missing
./configure
make
sudo make install
```

