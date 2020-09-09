
    This program generates hybrid coordinate interface profiles
    from a few given parameters. The hybrid coordinate is intended to start
    out at the bottom in pure sigma and end up at the top in pure pressure,
    with a smooth transition in between. The pressure thickness is close to
    quadratic in pressure, with maximum thicknesses in the middle of the domain.
    The coordinate pressure will have continuous second derivatives in level.
 
    The hybrid coordinate is returned in terms of vectors AK and BK, where
    the interface pressure is defined as A+B*ps where ps is surface pressure
    Thus A=0 in regions of pure sigma and B=0 in regions of pure sigma.
    At the bottom, A(0)=0 and B(0)=1 so that surface pressure is the bottom
    boundary condition, while at the top, A(levs)=ptop and B(levs)=0 so that
    the constant top pressure (which can be zero) is the top boundary condition.
 
    Input argument list:
      levs     integer number of levels 
      lupp     integer number of levels below pupp
      pbot     real nominal surface pressure (Pa) 
      psig     real nominal pressure where coordinate changes
               from pure sigma (Pa) 
      ppre     real nominal pressure where coordinate changes
               to pure pressure (Pa) 
      pupp     real nominal pressure where coordinate changes
               to upper atmospheric profile (Pa)
      ptop     real pressure at top (Pa)
      dpbot    real coordinate thickness at bottom (Pa) 
      dpsig    real thickness of zone within which coordinate changes
               to pure sigma (Pa) 
      dppre    real thickness of zone within which coordinate changes
               to pure pressure (Pa) 
      dpupp    real coordinate thickness at pupp (Pa) 
      dptop    real coordinate thickness at top (Pa) 
 
    Outputs a text file containing the 'ak' and 'bk' values -
    bottom to top.  The conversion to pressure is:

      pressure = ak + (bk * pbot)
 
    Graphical description of parameters and zones:
      ptop  ---  -----  ----------------------
            ...  dptop
            ---         zone U (upper atmos)
            ...
      pupp  ---  -----  ----------------------
            ...  dpupp
            ---  -----
            ...         zone P (pure pressure)
            ---
            ...
      ppre  ---  -----  ----------------------
            ...
            ---  dppre  zone T1 (transition 1)
            ...
            ---  -----  ----------------------
            ...
            ---
            ...         zone T2 (transition 2)
            ---
            ...
            ---  -----  ----------------------
            ...
            ---  dpsig  zone T3 (transition 3)
            ...
      psig  ---  -----  ----------------------
            ...
            ---  -----  zone S (pure sigma)
            ...  dpbot
      pbot  ---  -----  ----------------------
 
    For details on the procedure, see the program prolog.
