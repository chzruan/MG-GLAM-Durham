!  Stand-in for RanluxLevel() when ranlux_lags.f90 is built against a
!  PMP2mod_random.f90 from before 2026-10-19 (make ... ORIG=1): that module
!  has no such function.  Returns 2, the level every caller asked for then;
!  the old module ignored it and drew at level 0 whatever the value.
      Integer*4 Function RanluxLevel()
      RanluxLevel = 2
      End Function RanluxLevel
