import copy
import numpy as np
from prettytable import PrettyTable


class SimplexMethod:

    # в функции __init__ создаем матрицу из предоставленных данных, которая будет символизировать симплекс таблицу и с
    # через нее же будем находить опорную симплекс таблицу и оптимальные решения

    def __init__(self, A, c, b):
        temp = A
        for i in range(len(temp)):
            temp[i].append(b[i])
        temp.append(c)
        temp[len(temp) - 1].append(0)

        self.SimplexTable = np.matrix(temp)

        self.RowNum, self.ColNum = self.SimplexTable.shape
        self.SimplexTable[self.RowNum - 1] *= -1

        self.OldSimplexTable = copy.deepcopy(self.SimplexTable)

        self.IdxRow = 0
        self.IdxCol = 0

        self.Headers = [''] + [f'x{i}' for i in range(1, self.ColNum)] + ['C']
        self.BaseCol = [[f'u{i}'] for i in range(1, self.RowNum)] + [['F']]

    def get_pivot(self):
        return self.SimplexTable[self.IdxRow, self.IdxCol]

    def jordan_transform(self):

        pivot = self.get_pivot()

        self.SimplexTable[self.IdxRow] /= pivot

        self.SimplexTable[:, self.IdxCol] /= -pivot

        self.SimplexTable[self.IdxRow, self.IdxCol] *= -1

        for i in range(self.RowNum):
            for j in range(self.ColNum):

                if i == self.IdxRow or j == self.IdxCol:
                    continue
                else:

                    self.SimplexTable[i, j] = (pivot * self.SimplexTable[i, j] -
                                               self.OldSimplexTable[self.IdxRow, j] *
                                               self.OldSimplexTable[i, self.IdxCol]) / pivot
        self.OldSimplexTable = copy.deepcopy(self.SimplexTable)

    def check_free_members(self):
        flag = True
        for i in range(self.RowNum - 1):
            if self.SimplexTable[i, self.ColNum - 1] < 0:
                flag = False
                break
        return flag

    def ref_perm_col(self):
        if self.SimplexTable[0:self.RowNum - 1, self.ColNum - 1].min() < 0:
            self.IdxCol = self.SimplexTable[0:self.RowNum - 1, self.ColNum - 1].argmin()

    def perm_row(self):
        min_row_arg = None

        for i in range(self.RowNum - 1):
            flag = self.SimplexTable[i, self.IdxCol] != 0
            if min_row_arg is None and flag and self.SimplexTable[i, self.ColNum - 1] / self.SimplexTable[
                i, self.IdxCol] > 0:
                min_row_arg = self.SimplexTable[i, self.ColNum - 1] / self.SimplexTable[i, self.IdxCol]
                self.IdxRow = i
            elif flag and 0 < self.SimplexTable[i, self.ColNum - 1] / self.SimplexTable[i, self.IdxCol] <= min_row_arg:
                min_row_arg = self.SimplexTable[i, self.ColNum - 1] / self.SimplexTable[i, self.IdxCol]
                self.IdxRow = i

    def swap(self):
        temp = self.BaseCol[self.IdxRow]
        self.BaseCol[self.IdxRow] = [self.Headers[self.IdxCol + 1]]
        self.Headers[self.IdxCol + 1] = temp[0]

    def ref_solution(self):
        iteration = 1
        self.print_table()

        while not self.check_free_members():
            self.ref_perm_col()
            self.perm_row()
            self.jordan_transform()
            print()

            self.swap()

            self.print_table()
            print(f'iteration: {iteration}')
            print("permission col: ", self.IdxCol)
            print("permission row: ", self.IdxRow)
            iteration += 1

        self.search_optimal()

    def check_obj_func(self):
        flag = True
        for i in range(self.ColNum - 1):
            if self.SimplexTable[self.RowNum - 1, i] < 0:
                flag = False
                break
        return flag

    def perm_col(self):
        if self.SimplexTable[self.RowNum - 1, 0:self.ColNum - 1].min() < 0:
            self.IdxCol = self.SimplexTable[self.RowNum - 1, 0:self.ColNum - 1].argmin()

    def search_optimal(self):
        iteration = 1
        self.print_table()

        while not self.check_obj_func():
            self.perm_col()
            self.perm_row()
            self.jordan_transform()
            print()

            self.swap()

            self.print_table()
            print(f'iteration: {iteration}')
            print("permission col: ", self.IdxCol)
            print("permission row: ", self.IdxRow)
            iteration += 1
        print()
        print("answer F = ", round(self.SimplexTable[-1, -1], 2))

    def print_table(self):
        repr_table = PrettyTable()
        repr_table.field_names = self.Headers
        for i in range(self.RowNum):
            if i == self.RowNum - 1:
                repr_table.add_row(self.BaseCol[len(self.BaseCol) - 1] + np.array(self.SimplexTable)[i].tolist())
            else:
                repr_table.add_row(self.BaseCol[i] + np.array(self.SimplexTable)[i].tolist())
        return print(repr_table)


