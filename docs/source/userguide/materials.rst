.. _materials:

---------
Materials
---------

Material is just a simple compoistion of nuclides with concentration will be
changed over the time of the simulation process. To ensure right solution
some attributes should be added to material like ``volume``[cm^3], ``power`` if
materials is on induced radiation. The :ref:`main block <material_xml>` 
constitutes with two arrayes one with name of nuclides and another with nuclear
conctrations multiplied by 10^-24. To uncertainties analysis user can provide
initial uncertainties for every nuclide in according array ``<dconc>``.  

