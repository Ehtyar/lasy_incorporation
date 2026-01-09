(TeX-add-style-hook "beamerouterthemetud"
 (lambda ()
    (TeX-add-symbols
     '("datecity" 1)
     '("professur" 1)
     '("institut" 1)
     '("fachrichtung" 1)
     '("einrichtung" 1)
     "tudbeamer"
     "insertdatecity"
     "inserttotalpagenumber"
     "setbeamertemplates"
     "tudbeamerslidename"
     "tudbeamerofname"
     "frametitle"
     "noexpand")
    (TeX-run-style-hooks
     "calc")))

