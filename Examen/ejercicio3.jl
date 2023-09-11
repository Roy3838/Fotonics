

# Se tiene un material descrito por 
E = (2x + y + 2z) (E_0/sqrt(2))

# Y se tiene una susceptibilidad descrita por 
e_r = [1.5 0.5 0; 0.5 1.5 0; 0 0 1]


P = e_0*e_r*E
# Se quiere saber el angulo walk off, o el angulo entre D y E 
D = E_0*e_0*[4;2;2]

# Tambien se quiere saber los indices de refraccion dados en el sistema de ejes principal, que son los eigenvectores de e_r
