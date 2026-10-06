# Figure style

`fl_style.py` is the project's figure style module, copied from the FastaLake code
repository (`figure_style/fl_style.py`, sha256
`972011c2219b49e8a61b521f327863f086e398786cdaff1f9903c1c93d4729e5`). One block was
changed: the canonical copy imported two lab styleguide modules from a directory named by
`FL_STYLE_DIR`. The three symbols it used are reproduced in `_house.py` with the same
behaviour, so nothing outside this repository is needed except the typeface.

Usage in a notebook or `panels.py`:

    from style import fl_style as fl
    fl.init()

## Fonts

The figures are set in Söhne (Klim Type Foundry). The font is licensed and is not in this
repository. To reproduce the published panels pixel for pixel, place `soehne-buch.otf`,
`soehne-kraftig.otf` and `soehne-halbfett.otf` in a directory and set

    export FASTALAKE_FONT_DIR=/path/to/that/directory

(`FL_STYLE_DIR` is accepted as an alias.) Without it, `fl_style` warns once and text is
drawn in the matplotlib default face (DejaVu Sans). Colours, sizes, line weights and
layout are unchanged, but text pixels differ, so pixel checks of panels with text are only
meaningful with the fonts present. Continuous integration has no fonts and does not run
pixel checks.
