# FilterString examples

mzML files converted by msconvert and ThermoRawFileParser include the Thermo filter string that hold great information.
It looks like this:

<cvParam cvRef="MS" accession="MS:1000512" value="FTMS + c NSI d Full ms2 697.37@hcd32.00 [180.00-1405.00]" name="filter string" />

Here are some examples:

- FTMS + p NSI Full ms [350.00-1500.00]
High mass accuracy precursor: ms and FTMS

- FTMS + c NSI d Full ms2 697.37@hcd32.00 [180.00-1405.00]
High mass accuracy (FTMS) HCD = 

- ITMS + c NSI r d sa Full ms2 470.93@etd33.33 [50.00-1425.00]
Low mass accuracy fragmentation (ITMS) with ETD

- ITMS + c NSI r d Full ms2 509.2148@hcd30.00 [110.0000-1029.0000]
Low mass accuracy (ITMS) HCD here apparently

- FTMS + c NSI d sa Full ms2 557.2769@etd116.52@cid30.00 [120.0000-1125.0000]
High mass accuracy (FTMS) ETD + CID

- ITMS + c NSI d w Full ms2 857.62@cid30.00 [225.00-870.00]
Plain low mass accuracy (ITMS) CID

- At PeptideAtlas, we track the following types and filter strings:
1   HR IT CID       LTQ FT instruments: FTMS + c NSI d Full ms2 697.37@cid32.00 [180.00-1405.00]
2   HR IT ETD       LTQ FT instruments: FTMS + c NSI d Full ms2 697.37@etd32.00 [180.00-1405.00]
3   HR HCD          FTMS + c NSI d Full ms2 697.37@hcd32.00 [180.00-1405.00]
4   HR Q-TOF         (not Thermo; no filter string)
5   LR IT CID       ITMS + c NSI d w Full ms2 857.62@cid30.00 [225.00-870.00]
6   LR IT ETD       ITMS + c NSI r d sa Full ms2 470.93@etd33.33 [50.00-1425.00]
7   LR HCD          ITMS + c NSI r d Full ms2 509.2148@hcd30.00 [110.0000-1029.0000]                * (HCD IT)
8   HR ETciD        FTMS + c NSI d sa Full ms2 557.2769@etd116.52@cid30.00 [120.0000-1125.0000]     Have we observed this? Is it possible?
9   LR ETciD        ITMS + c NSI r d sa Full ms2 984.7690@etd95.59@hcd25.00 [110.0000-1980.0000]    * (EThcD_IT) (ETcaD_IT)
10  HR EThcD        FTMS + p NSI d sa Full ms2 982.3831@etd95.59@hcd25.00 [110.0000-1975.0000]      * (EThcD OT)
