using MyLib;

const double tolerance = 0.001;
const int maxIterations = 10000;

double F(Matrix point)
{
    double x1 = point[0, 0];
    double x2 = point[0, 1];
    return x1 * x1 + x2 * x2 - 4 * x1 - 6 * x2;
}

Matrix Grad(Matrix point)
{
    double x1 = point[0, 0];
    double x2 = point[0, 1];
    return new Matrix(new double[1, 2] { { 2 * x1 - 4, 2 * x2 - 6 } });
}

// Проверяем находится ли точка в X
bool Check(Matrix point)
{
    double x1 = point[0, 0];
    double x2 = point[0, 1];
    if (Math.Abs(2 * x1 + x2 - 2.0) < tolerance) return true;
    return false;
}

// Находим проекцию точки на прямую 2 * x1 + x2 = 2
Matrix GetProj(Matrix point)
{
    double x1 = point[0, 0];
    double x2 = point[0, 1];

    double b = x2 - 0.5 * x1;
    double res1 = (2 - b) / 2.5;
    double res2 = 2 - 2 * res1;
    
    return new Matrix(new double[1, 2] { { res1, res2 } });
}

// Решение второй задачи методом проекции градиента
Matrix Solve(
    Matrix start
)
{
    Matrix x = start;
    Matrix xPast = x;
    double alpha = 1;
    for (int i = 0; i < maxIterations; i++)
    {
        // Подбираем alpha
        for (int j = 0; j < maxIterations; j++)
        {
            x = GetProj(x - alpha * Grad(x));
            
            // Условие монотонности
            if (F(x) < F(xPast))
            {
                break;
            }
            x = xPast;
            alpha /= 2;
        }
        
        // Выход из цикла
        if (Methods.EuclidNorm(x - xPast) < tolerance)
        {
            break;
        }

        xPast = x;
    }

    return x;
}

// Начальная точка (1, 0)
Matrix st = new Matrix(new double[1, 2] { { 1, 0 } });

Matrix result = Solve(st);
result.Print();
Console.WriteLine($"Минимальное значение: {F(result)}");
