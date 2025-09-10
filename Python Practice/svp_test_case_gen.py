from fpylll import *
from itertools import product

bits = [8, 16, 32]
#types = ["uniform", "ntrulike", "ntrulike2"]
types = ["uniform"]
dims_x_counts = [
    #(0,     1       ),
    #(1,     1_000   ),
    (2,     1_000   ),
    (3,     1_000   ),
    (4,     100     ),
    (5,     100     ),
    (6,     10      ),
    (7,     5       ),
    (8,     3       ),
    (9,     2       ),
    (10,    1       ),
]

for b, (d, c), t in product(bits, dims_x_counts, types):
    for i in range(c):
        mat = IntegerMatrix(d, d)
        mat.randomize(t, bits=b)
        name = f"{d}-{i}_{b}_{t}"

        with open(f"examples/gen/{name}.svp", "w+") as f:
            mat_str = str([[c for c in l] for l in mat])
            f.write(mat_str)
        SVP.shortest_vector(mat)
        with open(f"examples/gen/{name}.svs", "w+") as f:
            f.write(str(mat[0].norm()))