

function critical_lambda(d,n1,n2,n3)
    AN = sqrt(n1*n1 - n2*n2)
    A = 2 * d * AN
    B = tand((n1*n1*(1/AN) - n3*n3 )/(n1*AN))/pi
    return A / B
end

# critical_lambda for d = 0.5 um ns are sellmeier

function get_n(A,B,C, lambda)
    t1 = A[0] +
    (B[0]*lambda*lambda)/(lambda*lambda - C[0]) +
    (B[1]*lambda*lambda)/(lambda*lambda - C[1]) +
    (B[2]*lambda*lambda)/(lambda*lambda - C[2])
    return sqrt(t1 + t2 + t3 + t4)
end

# Si3N4
A= [1]
B= [3.0249, 40314, 0]
C= [0.018311707 , 1.5372081e6, 0]

# SiO2
A = [1]
B = [0.6961663, 0.4079426, 0.8974794]
C = [0.004679148, 0.01351206, 97.934002]

# Silicio
A = [1]
B = [10.861286, -1.6232027e6, 0]
C = [0.1273519, 8.1569e23, 0]





