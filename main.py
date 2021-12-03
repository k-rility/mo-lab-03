import from_input
import SimplexMethod
import copy

if __name__ == '__main__':
    A, c, b = from_input.get_data(open("input.txt", "r"))
    res = SimplexMethod.SimplexMethod(A, c, b)
    res.ref_solution()
