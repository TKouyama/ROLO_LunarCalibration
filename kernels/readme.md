In this directory, programs expect below directories:

- lsk/
- spk/planets
- fk/planets/
- fk/satellites/
- pck/

As in luna_pos.tm, below kernels are used in programs:

      \begindata
      KERNELS_TO_LOAD = ( 'lsk\naif0012.tls.pc',
                          'spk\planets\de430.bsp',
                          'fk\planets\earth_assoc_itrf93.tf',
                          'fk\planets\RSSD0002.TF',
                          'fk\satellites\moon_080317.tf',
                          'pck\moon_pa_de421_1900-2050.bpc',
                          'pck\earth_720101_070426.bpc',
                          'pck\earth_000101_160620_160330.bpc',
                          'pck\earth_070425_370426_predict.bpc',
                          'pck\pck00010.tpc'
                         )
