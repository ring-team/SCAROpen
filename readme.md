# SCAROpen
## Sealed model Creation and Automatic Repair
## By RING


SCAROpen is a [RINGMesh][ringmesh]-dependant library (shipped with some executables) derived from [SCAR](https://www.ring-team.org/technologies/181-scar), a code mainly developed by Pierre Anquez during his PhD (2016-2019), within the [RING-team][ring], GeoRessources, University of Lorraine, France. His PhD was financed by the industrial and academic sponsors of the [RING Consortium][consortium] managed by ASGA.

The RING Consortium Steering Committee is deeply thanked for allowing the creation of this open-source repository. 

## License

Copyright (c) 2021 ASGA and Universite de Lorraine. All Rights Reserved.

This program (SCAROpen) is the open version of the SCAR program developed 
in the frame of the RING project managed by ASGA and Universite de Lorraine. 

It is distributed under a dual licensing scheme:

1. RING Consortium License
Members of the RING-COCAD Consortium may only use this file in
accordance with the terms of described in the GOCAD Advancement Agreement, 
without the prior written authorization of the ASGA.
Licencee agrees to attach or embed this Notice on all copies 
of the program, including partial copies or modified versions thereof.
Please use: contact at ring dash team dot org, for more information. 

2. GNU General Public License Usage
Alternatively, this file may be used under the terms of the 
GNU General Public license version 3. The licenses are as published by 
the Free Software Foundation and appearing in the file 
https://www.gnu.org/licenses/gpl-3.0.html

## Contact

RING team: contact at ring dash team dot org

## Features

Repair or simplify a 2D GeoModel (cross-sections, maps) as described in [Anquez et al. 2019][papierCRG] and [Anquez 2019][these] (PHD thesis in French).

## Articles and references

This reposity is dedicated to public the algorithms used in:

Anquez, P., Glinsky, N., Cupillard, P., & Caumon, G. (2022). Impacts of geometric model simplifications on wave propagation - Application to ground motion simulation in the lower Var valley basin (France), Geophysical Journal International, 229, 110–137.
[Click here][papierGJI]

Another article and a PhD thesis (in French) present more in details the algorithms:

Anquez, P., Pellerin, J., Irakarama, M., Cupillard, P., Lévy, B., & Caumon, G. (2019). Automatic correction and simplification of geological maps and cross-sections for numerical simulations. Comptes Rendus Geoscience, 351(1), 48-58.
[Click here][papierCRG]

Anquez, P. (2019). Correction et simplification de modèles géologiques par frontières: impact sur le maillage et la simulation numérique en sismologie et hydrodynamique (Doctoral dissertation, Université de Lorraine). [Click here][these]


## Disclamer

SCAROpen is derived from SCAR, a research code written by Pierre Anquez between 2016 and 2019. This code is research code which includes some unfinished or sub-optimal implementation, and probably several bugs. There has been no more development on this code in the RING team since 2019 at the moment we write theses lines (October 2021). Please do not hesitate to contact RING for any comment, question or information. 


[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [ringmesh]: <https://github.com/ringmesh/RINGMesh>
   [ring]: <https://www.ring-team.org/>
   [consortium]: <https://www.ring-team.org/consortium>
   [papierCRG]: <https://www.sciencedirect.com/science/article/pii/S1631071318301706>
   [these]: <https://hal.univ-lorraine.fr/tel-02330956/document>
   [papierGJI]: <https://doi.org/10.1093/gji/ggab447>
