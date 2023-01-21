using Plots

# CONSTANT VALUES

# constant G
const G = 6.67408e-11
# domain of the problem
const range = [0, 3]
# input of parameter for the number of divisions !!!
print("Enter the number of divisions: ")
const n = parse(Int64, readline()) # input 
# distance between each division
const h = (range[2] - range[1])/n
# weights 1 and 2 for the Gauss Legendre Quadrature (GLQ) - 4 points formula
const w1 = (18 + sqrt(30))/36
# weights 3 and 4 for the Gauss Legendre Quadrature (GLQ) - 4 points formula
const w2 = (18 - sqrt(30))/36
# points 1 and 2 for the Gauss Legendre Quadrature (GLQ) - 4 points formula
const x12 = sqrt((3/7)-(2/7)*sqrt(6/5))
# point 3 for the Gauss Legendre Quadrature (GLQ) - 4 points formula
const x3 = sqrt((3/7)+(2/7)*sqrt(6/5))
# point 4 for the Gauss Legendre Quadrature (GLQ) - 4 points formula
const x4 = -sqrt((3/7)+(2/7)*sqrt(6/5))

# function rho based on the given data in the problem
function rho(x::Number)
    if x > 1 && x < 2
        return 1
    else
        return 0
    end
end

# Auxiliary function to calculate the new x for the Gauss Legendre Quadrature (GLQ) - 4 points formula
function calculateNewX(x::Number, r1::Number, r2::Number)
    return (r2 - r1)/2 * x + (r1 + r2)/2
end


# Numerical Integration Gauss Legendre Quadrature (GLQ) - 4 points formula
function integral(f, r1::Number, r2::Number) # f is the function, r1 and r2 are the limits of integration
    # return the result of the integration
    return (r2 - r1)/2 * (w1*f(calculateNewX(x12, r1, r2)) + w1*f(calculateNewX(x12, r1, r2)) + w2*f(calculateNewX(x3, r1, r2)) + w2*f(calculateNewX(x4, r1, r2)))
end

function e_i(i::Int64)
    return function(x)
        if x > h*i && x < h*(i+1)
            return (h*(i+1) - x)/h
        elseif x > h*(i-1) && x <= h*i
            return (x - h*(i-1))/h
        else
            return 0
        end
    end
end

function e_i_prime(i::Int64)
    return function(x)
        if x > h*i && x < h*(i+1)
            return -1/h
        elseif x > h*(i-1) && x <= h*i
            return 1/h
        else
            return 0
        end
    end
end

function FEM(n) 
    # Preparing matrices
    A = zeros(n-1, n-1)
    b = zeros(n-1)

    # Filling the matrix A

    # Filling diagonal
    f1(x) = e_i_prime(1)(x) * e_i_prime(1)(x)
    tmp = -integral(f1, h*0, h*2)
    for i in 1:n-1
        A[i, i] = -tmp
    end
    # Filling symetrical triangles
    for i in 1:n-2
        f2(x) = e_i_prime(i)(x)*e_i_prime(i+1)(x)
        tmp = -integral(f2, h*i, h*(i+1))
        A[i, i+1] = tmp
        A[i+1, i] = tmp
    end

    # Filling the vector b
    for i in 1:n-1
        f3(x) = e_i_prime(i)(x)*rho(x)
        f4(x) = -(1/3) * e_i_prime(i)(x)
        b[i] = 4*pi*G*(integral(f3, h*(i-1), h*i) + integral(f3, h*i, h*(i+1))) + integral(f4, h*(i-1), h*(i)) + integral(f4, h*i, h*(i+1))
    end

    # Solving the system
    x = A\b

    # Plotting the solution

    # Returning the solution
    return function(x)
        solution = 5 - (x/3)
        for i = 1:n-1
            solution += B[i]*e_i(i)(x)
        end
        return solution  
    end
end

# Running the code
FEM(n)
# f(x) = e_i_prime(1)(1.5)
# f(2)