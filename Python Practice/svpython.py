# Shortest Vector Problem for 10 dimensions or less using LLL algorithm and Kannan's algorithm
import sys 

def kahan(list):
    summation = 0
    error = 0
    for i in list:
        y = i - error
        temp = summation + y
        error = (temp - summation) - y
        summation = temp
    return summation

def dotproduct(v1, v2):
    print("Vector 1: ", v1)
    print("Vector 2: ", v2)
    # print("Dot product: ", sum(x * y for x, y in zip(v1, v2)))
    # return sum(x * y for x, y in zip(v1, v2))
    kahanana = kahan([v1[i]*v2[i] for i in range(len(v1))])
    print("Dot product with Kahan: ", kahanana)
    return kahanana

def subvec(v1, v2):
    result = []
    for i in range(len(v1)):
        result.append(v1[i] - v2[i])
    print("Subvec: ", result)
    return result

def multvect(v1, scalar):
    print("Multiply vector!")
    if v1 is None or not isinstance(v1, list):
        return []
    print("V1 is a non-empty list")
    print("V1 is: ", v1)
    result = []
    print("V1 length: ", len(v1))
    for i in range(len(v1)):
        print("Result now: ", result)
        print("V1[i]: ", v1[i])
        print("Scalar: ", scalar)
        print("V1[i] * Scalar: ", (v1[i] * scalar))
        result.append(v1[i] * scalar)
    print("Multvect: ", result)
    return result

def check_linearly_independent(basis):
    matrix = [list(vec) for vec in basis]
    n = len(matrix)
    for i in range(n):
        if matrix[i][i] == 0:
            for j in range(i+1, n):
                if matrix[j][i] != 0:
                    matrix[i], matrix[j] = matrix[j], matrix[i]
                    break
            else:
                return False
        for j in range(i+1, n):
            ratio = matrix[j][i] / matrix[i][i]
            matrix[j] = [matrix[j][k] - ratio * matrix[i][k] for k in range(n)]
    return True

def gs(basis, m): # Gram Schmidt orthogonalization
    for i in range(m):
        for j in range(i):
            print("Gram Schmidt l41")
            dpjj = dotproduct(basis[j], basis[j])
            if dpjj == 0:
                continue
            print("Gram Schmidt l45")
            mu = round(dotproduct(basis[i], basis[j]) / dpjj, 100)
            basis[i] = subvec(basis[i], multvect(basis[j], mu))
    return basis

def sr(basis, m): # Size reduction
    for i in range(m):
        for j in range(i):
            print("Size reduction l56")
            dpjj = dotproduct(basis[j], basis[j])
            if dpjj == 0:
                continue
            print("Size reduction l60")
            mu = round(dotproduct(basis[i], basis[j]) / dpjj, 100)
            basis[i] = subvec(basis[i], multvect(basis[j], mu))
    return basis

def lll_reduction(basis):
    blocksize = 2
    reduction = 0.75
    n = len(basis[0])
    m = len(basis)

    if not check_linearly_independent(basis):
        return basis

    for i in range(1, m):
        for j in range(i-1, -1, -1):
            # Gram-Schmidt orthogonalization
            print("LLL l73")
            dpjj = dotproduct(basis[j], basis[j])
            if dpjj != 0:
                print("LLL l76")
                mu = round(dotproduct(basis[i], basis[j]) / dpjj, 100)
                basis[i] = subvec(basis[i], multvect(basis[j], mu))   
        # Size reduction with Lovasz condition
        sqnorm = norm_sq(basis[i])
        print("LLL l81")
        dpim1 = dotproduct(basis[i-1], basis[i-1])
        if i > 1 and dpim1 != 0 and abs(dotproduct(basis[i], basis[i-1]) / dpim1) > reduction / 4 and sqnorm >= reduction * dpim1:
            print("LLL l84")
            mu = round(dotproduct(basis[i], basis[i-1]) / dpim1, 100)
            basis[i] = subvec(basis[i], multvect(basis[i-1], mu))
            if norm_sq(basis[i]) == 0:
                continue
            basis[i], basis[i-1] = basis[i-1], basis[i]
            i -= 1
        basis = gs(basis, len(basis))
    return basis

def copy(list):
    newlist = []
    for i in range(0, len(list)):
        newlist.append(list[i])
    return newlist

def candidate_gen(basis):
    # Generating candidates using BKZ for shortest vector
    newb = copy(basis)
    n = len(basis[0]) # We're taking the length of the first vector since they should all be the same dimensionality
    m = len(basis) # Number of vectors in lattice basis
    blocksize = 2 # I used 2 before so I'll use 2 again and change if necessary

    while blocksize < m:
        basis = gs(basis, m)
        basis = sr(basis, m)

        for i in range(m):
            for j in range(i+1, min(i+blocksize, m)):
                print("Candidate gen l112 (Checking dpii)")
                dpii = dotproduct(basis[i], basis[i])
                if dpii == 0:
                    continue
                print("Candidate gen l116")
                mu = round(dotproduct(basis[j], basis[i]) / dpii, 100)
                print("CG mu: ", mu)
                basis[j] = subvec(basis[j], multvect(basis[i], mu)) 
                print("New CG Basis: ", basis)

        blocksize = min(2 * blocksize, m)
    
    return basis

def norm_sq(vector):
    result = 0
    if isinstance(vector, int):
        print("Vector was int. Norm squared: ", vector ** 2)
        return vector ** 2
    else:
        for i in range(0, len(vector)):
            result += vector[i] * vector[i]
        print("Norm squared: ", result)
        return result

def enum_neighbours(v, d=1):
    neighbours = []
    n = len(v)

    def enum_gen(i, vec, n):
        if i == n:
            neighbours.append(vec)
            return
        enum_gen(i+1, vec[:i] + [vec[i] + d] + vec[i+1:], n)
        enum_gen(i+1, vec[:i] + [vec[i] - d] + vec[i+1:], n)

    enum_gen(0, v, n)
    return neighbours

def kannan(candidates):
    print("Into Kannan now")
    print("Kannandidates: ", candidates)
    best_v = None
    best_sqnorm = float('inf')

    for candidate in candidates:
        current_v = candidate
        current_sqnorm = norm_sq(current_v)
        print("Current square norm: ", current_sqnorm)

        for neighbour in enum_neighbours(current_v):
            neighbour_sqnorm = norm_sq(neighbour)
            if neighbour_sqnorm < best_sqnorm:
                best_v = neighbour
                best_sqnorm = neighbour_sqnorm
                print("Best square norm: ", best_sqnorm)

    return best_sqnorm

def solve_svp(basis):
    if check_validity(basis) == False:
        return "Invalid basis"
    reduced_basis = lll_reduction(basis) # Reduce our basis using LLL
    candidates = candidate_gen(reduced_basis) # Generate candidates for shortest vector
    shortest_vector = kannan(candidates) # Pick the shortest vector from the candidates
    return shortest_vector

def check_validity(basis):
    if len(basis) == 0: # We can't have a 0-dimensional basis!
        print("Basis is empty")
        return False
    for i in range(0, len(basis)):
        if len(basis[i]) != len(basis): # If the length of each basis vector isn't the same as the number of basis vectors...
            print("Length of each vector basis is not the same as the number of basis vectors")
            return False                # ...then the basis is invalid.
    
    if check_linearly_independent(basis) == False:
        print("Basis is not linearly independent")
        return False
    
    return True # Otherwise we're all good. 
 
def build_basis(argv):
    init = []
    basis = []
    sublist = []
    for i in range(1, len(argv)): # The 0th argument is the program name.
        init.append(argv[i])
    
    for j in range(0, len(init)):
        init[j] = init[j].replace(',','')
        if '[' in init[j]:
            init[j] = init[j].replace('[','')
            sublist.append(int(init[j]))
        elif ']' in init[j]:
            init[j] = init[j].replace(']','')
            sublist.append(int(init[j]))
            basis.append(sublist)
            sublist = []
        else:
            sublist.append(int(init[j]))
    
    return basis

basis = build_basis(sys.argv)
print(solve_svp(basis))
