from fractions import Fraction

def get_matrix_input():
    n = int(input("정수 n 입력: "))
    print(f"{n}x{n} 행렬의 원소를 행 단위로 공백으로 구분하여 입력하세요:")
    matrix = []
    for _ in range(n):
        row = [Fraction(x) for x in input().split()]
        matrix.append(row)
    return n, matrix

def get_minor(m, i, j):
    return [row[:j] + row[j+1:] for row in (m[:i] + m[i+1:])]

def determinant(m):
    if len(m) == 1:
        return m[0][0]
    if len(m) == 2:
        return m[0][0]*m[1][1] - m[0][1]*m[1][0]
    det = Fraction(0)
    for c in range(len(m)):
        det += ((-1)**c) * m[0][c] * determinant(get_minor(m, 0, c))
    return det

def inverse_determinant(m):
    det = determinant(m)
    if det == 0:
        return None
    if len(m) == 1:
        return [[Fraction(1) / det]]
    cofactors = []
    for r in range(len(m)):
        cofactor_row = []
        for c in range(len(m)):
            minor = get_minor(m, r, c)
            cofactor_row.append(((-1)**(r+c)) * determinant(minor))
        cofactors.append(cofactor_row)
    adjugate = list(map(list, zip(*cofactors)))
    inverse = []
    for r in range(len(adjugate)):
        inverse_row = []
        for c in range(len(adjugate)):
            inverse_row.append(adjugate[r][c] / det)
        inverse.append(inverse_row)
    return inverse

def inverse_gauss_jordan(m):
    n = len(m)
    aug = [row[:] + [Fraction(1) if i == j else Fraction(0) for j in range(n)] for i, row in enumerate(m)]
    for i in range(n):
        if aug[i][i] == 0:
            for k in range(i+1, n):
                if aug[k][i] != 0:
                    aug[i], aug[k] = aug[k], aug[i]
                    break
            else:
                return None
        pivot = aug[i][i]
        for j in range(i, 2*n):
            aug[i][j] /= pivot
        for k in range(n):
            if k != i:
                factor = aug[k][i]
                for j in range(i, 2*n):
                    aug[k][j] -= factor * aug[i][j]
    inverse = [row[n:] for row in aug]
    return inverse

def compare_matrices(m1, m2):
    if m1 is None or m2 is None:
        return m1 == m2
    for r in range(len(m1)):
        for c in range(len(m1[0])):
            if m1[r][c] != m2[r][c]:
                return False
    return True

def print_matrix(m):
    for row in m:
        print(" ".join([str(x) for x in row]))

def main():
    n, matrix = get_matrix_input()
    
    det = determinant(matrix)
    if det == 0:
        print("\n오류: 역행렬이 존재하지 않습니다 (행렬식이 0입니다).")
        return
        
    inv_det = inverse_determinant(matrix)
    print("\n[방법 1] 행렬식을 이용한 역행렬 (분수 표현):")
    print_matrix(inv_det)
    
    inv_gj = inverse_gauss_jordan(matrix)
    print("\n[방법 2] 가우스-조던 소거법을 이용한 역행렬 (분수 표현):")
    print_matrix(inv_gj)
    
    print("\n[결과 비교]")
    if compare_matrices(inv_det, inv_gj):
        print("두 방법으로 계산한 역행렬이 동일합니다.")
    else:
        print("두 방법의 결과가 다릅니다.")

if __name__ == '__main__':
    main()