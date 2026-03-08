======================
Upcoming Release Notes
======================

---------
Bug Fixes
---------

- There was a bug where copying a DiffusionData instance in Python did not include the
  LeakageCorrections which may be present. This resulted in the copies not having a
  LeakageCorrections instance and gave incorrect results in nodal diffusion calculations.

