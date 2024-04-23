
include("ftdt.jl")
using .FTDT_LIB

#

########################################################################
# Se define la geometría en la ventana numérica principal
#######################################################################

# Se agregan objetos en el medio
casoMaterial = "vacio";
if casoMaterial=="vacio"
  idx = 1:je;
  idy = 1:ie;
  medioProp = epsz;

simular(medioProp, idx, idy)
