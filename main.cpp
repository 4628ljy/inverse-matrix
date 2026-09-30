#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// 행렬 출력
void printMatrix(const vector<vector<double>>& a)
{
    int n = a.size();

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double value = a[i][j];
            if (fabs(value) < 1e-10)
                value = 0;
            cout << value << "   ";
        }
        cout << '\n';
    }
}

// 전치행렬 구하기
vector<vector<double>> transposeMatrix(const vector<vector<double>>& a)
{
    int n = a.size();
    vector<vector<double>> result(n, vector<double>(n));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            result[j][i] = a[i][j];
        }
    }

    return result;
}

// 특정 행과 열을 제외한 소행렬 구하기
vector<vector<double>> getMinor(const vector<vector<double>>& a, int row, int col)
{
    int n = a.size();
    vector<vector<double>> minor;

    for (int i = 0; i < n; i++)
    {
        if (i == row) continue;

        vector<double> temp;

        for (int j = 0; j < n; j++)
        {
            if (j != col)
                temp.push_back(a[i][j]);
        }

        minor.push_back(temp);
    }

    return minor;
}

// 행렬식 계산
double determinant(const vector<vector<double>>& a)
{
    int n = a.size();

    if (n == 1)
    {
        return a[0][0];
    }

    if (n == 2)
    {
        return a[0][0] * a[1][1] - a[0][1] * a[1][0];
    }

    double result = 0;

    for (int i = 0; i < n; i++)
    {
        vector<vector<double>> minor = getMinor(a, 0, i);

        if (i % 2 == 0)
            result += a[0][i] * determinant(minor);
        else
            result -= a[0][i] * determinant(minor);
    }

    return result;
}

// 행렬식을 이용한 역행렬 계산
bool inverseByDeterminant(const vector<vector<double>>& a, vector<vector<double>>& inverse)
{
    int n = a.size();
    double det = determinant(a);

    if (fabs(det) < 1e-10)
    {
        return false;
    }

    inverse.assign(n, vector<double>(n));

    // 1차 행렬의 역행렬
    if (n == 1)
    {
        inverse[0][0] = 1.0 / a[0][0];
        return true;
    }

    // 2차 행렬의 역행렬
    if (n == 2)
    {
        inverse[0][0] = a[1][1] / det;
        inverse[0][1] = -a[0][1] / det;
        inverse[1][0] = -a[1][0] / det;
        inverse[1][1] = a[0][0] / det;

        return true;
    }

    // 3차, 4차 행렬의 여인수 행렬
    vector<vector<double>> cofactors(n, vector<double>(n));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            vector<vector<double>> minor = getMinor(a, i, j);

            double value = determinant(minor);

            if ((i + j) % 2 == 1)
                value = -value;

            cofactors[i][j] = value;
        }
    }

    // 여인수 행렬을 전치하여 수반행렬 구하기
    vector<vector<double>> adjugate = transposeMatrix(cofactors);

    // 행렬식으로 나누기
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            inverse[i][j] = adjugate[i][j] / det;
        }
    }

    return true;
}

// 가우스-조던 소거법을 이용한 역행렬 계산
bool inverseByGaussJordan(const vector<vector<double>>& a, vector<vector<double>>& inverse)
{
    int n = a.size();
    vector<vector<double>> temp(n, vector<double>(2 * n, 0));

    // 확장행렬 만들기
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            temp[i][j] = a[i][j];
        }
        temp[i][n + i] = 1;
    }

    for (int i = 0; i < n; i++)
    {
        // 현재 열에서 절댓값이 가장 큰 원소의 행 선택
        int pivotRow = i;

        for (int j = i + 1; j < n; j++)
        {
            if (fabs(temp[j][i]) > fabs(temp[pivotRow][i]))
            {
                pivotRow = j;
            }
        }

        // 피벗을 찾을 수 없으면 역행렬이 없음
        if (fabs(temp[pivotRow][i]) < 1e-10)
            return false;

        // 현재 행과 피벗 행 교환
        if (pivotRow != i)
        {
            vector<double> row = temp[i];
            temp[i] = temp[pivotRow];
            temp[pivotRow] = row;
        }

        // 대각 원소를 1로 만들기
        double pivot = temp[i][i];

        for (int j = 0; j < 2 * n; j++)
        {
            temp[i][j] /= pivot;
        }

        // 나머지 행의 현재 열 원소를 0으로 만들기
        for (int j = 0; j < n; j++) {
            if (j == i) continue;

            double value = temp[j][i];

            for (int k = 0; k < 2 * n; k++)
            {
                temp[j][k] -= value * temp[i][k];
            }
        }
    }

    // 확장행렬의 오른쪽 부분을 가져오기
    inverse.assign(n, vector<double>(n));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            inverse[i][j] = temp[i][n + j];
        }
    }

    return true;
}


// 두 역행렬 비교
bool compareMatrix(const vector<vector<double>>& a, const vector<vector<double>>& b)
{
    int n = a.size();

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            // 실수 계산 오차를 고려
            if (fabs(a[i][j] - b[i][j]) > 1e-10)
                return false;
        }
    }

    return true;
}


// 두 행렬의 곱 계산
vector<vector<double>> multiplyMatrix(const vector<vector<double>>& a, const vector<vector<double>>& b)
{
    int n = a.size();
    vector<vector<double>> result(n, vector<double>(n, 0));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n; k++)
            {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    return result;
}

// 단위행렬인지 확인
bool isIdentityMatrix(const vector<vector<double>>& a)
{
    int n = a.size();

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                // 대각 원소는 1인지 확인
                if (fabs(a[i][j] - 1) > 1e-10)
                    return false;
            }
            else
            {
                // 나머지 원소는 0인지 확인
                if (fabs(a[i][j]) > 1e-10)
                    return false;
            }
        }
    }

    return true;
}

int main()
{
    // 행렬입력
    int n;

    cout << "정방행렬의 차수를 입력하세요: ";
    cin >> n;
    cout << '\n';

    vector<vector<double>> a(n, vector<double>(n));

    for (int i = 0; i < n; i++)
    {
        cout << i + 1 << "행: ";
        for (int j = 0; j < n; j++)
        {
            cin >> a[i][j];
        }
    }
    vector<vector<double>> inverse1;
    vector<vector<double>> inverse2;

    bool result1 = inverseByDeterminant(a, inverse1);
    bool result2 = inverseByGaussJordan(a, inverse2);

    cout << "\n행렬식으로 구한 역행렬:\n";

    if (result1)
        printMatrix(inverse1);
    else
        cout << "역행렬이 존재하지 않습니다.\n";


    cout << "\n가우스-조던 소거법으로 구한 역행렬:\n";

    if (result2)
        printMatrix(inverse2);
    else
        cout << "역행렬이 존재하지 않습니다.\n";
    cout << '\n';

    if (result1 && result2)
    {
        if (compareMatrix(inverse1, inverse2))
            cout << "두 방법의 결과가 동일합니다.\n";
        else
            cout << "두 방법의 결과가 다릅니다.\n";
    }
    else
    {
        cout << "역행렬이 존재하지 않습니다.\n";
    }


    // 행렬식으로 구한 역행렬 검증
    if (result1)
    {
        vector<vector<double>> product1 = multiplyMatrix(a, inverse1);

        cout << "\n[행렬식 방법]검증\n";
        printMatrix(product1);

        if (isIdentityMatrix(product1))
            cout << "단위행렬과 일치합니다. 검증 성공!\n";
        else
            cout << "단위행렬과 일치하지 않습니다. 검증 실패!\n";
    }

    // 가우스-조던 소거법으로 구한 역행렬 검증
    if (result2)
    {
        vector<vector<double>> product2 = multiplyMatrix(a, inverse2);

        cout << "\n[가우스-조던 방법]검증\n";
        printMatrix(product2);

        if (isIdentityMatrix(product2))
            cout << "단위행렬과 일치합니다. 검증 성공!\n";
        else
            cout << "단위행렬과 일치하지 않습니다. 검증 실패!\n";
    }

    return 0;
}