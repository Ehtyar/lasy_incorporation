(TeX-add-style-hook "beamerstyle"
 (lambda ()
    (TeX-add-symbols
     "dinfamily")
    (TeX-run-style-hooks
     "beamerarticle"
     "babel"
     "fontenc"
     "T1"
     "inputenc"
     "utf8"
     "latex2e"
     "tudbook10"
     "tudbook"
     "ngerman"
     "oneside")))

