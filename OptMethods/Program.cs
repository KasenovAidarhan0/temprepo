using MethodsLib;

double[,] startPoint = { { 0, 0 } };
double epsilon = 0.001;
double alpha = 1;

double MyFunction(double[,] x)
{
    double x1 = x[0, 0];
    double x2 = x[0, 1];
    return 4 * x1 * x1 + x2 * x2 - x1 * x2 + 7 * x1 + 3 * x2;
}

double[,] Gradient(double[,] x)
{
    double x1 = x[0, 0];
    double x2 = x[0, 1];
    double[,] res = new double[1, 2];
    res[0, 0] = 8 * x1 - x2 + 7;
    res[0, 1] = 2 * x2 - x1 + 3;
    return res;
}

double[,] mR = { { 7, 3 } };
double[,] mQ =
{
    { 8, -1 },
    { -1, 2 }
};

double[,] NewGradient(double[,] x)
{
    double x1 = x[0, 0];
    double x2 = x[0, 1];
    double[,] res = new double[1, 2];
    res[0, 0] = 2 * x1 + x2 - 7;
    res[0, 1] = x1 + 4 * x2 - 7;
    return res;
}

double[,] newQ =
{
    { 2, 1 },
    { 1, 4 }
};

// Метод покоординатного спуска
double[,] min1 = Methods.CoordinateDescent(MyFunction, startPoint, epsilon, alpha);
// Метод Гаусса-Зейделя
double[,] min2 = Methods.GaussSeidel(mQ, mR, startPoint, epsilon);
// Метод градиентного спуска
double[,] min3 = Methods.GradientDescentDivide(MyFunction, Gradient, startPoint, epsilon);
// Метод наискорейшего спуска
double[,] min4 = Methods.GradientDescent(Gradient, mQ, startPoint, epsilon);
// Метод сопряженных градиентов
double[,] min5 = Methods.ConjugateDescent(NewGradient, newQ, startPoint, epsilon);

Methods.DisplayMatrixLine(min1, "Метод покоординатного спуска");
Console.WriteLine($"Минимальное значение: {MyFunction(min1)}");
Console.WriteLine();
Methods.DisplayMatrixLine(min2, "Метод Гаусса-Зейделя");
Console.WriteLine($"Минимальное значение: {MyFunction(min2)}");
Console.WriteLine();
Methods.DisplayMatrixLine(min3, "Метод градиентного спуска");
Console.WriteLine($"Минимальное значение: {MyFunction(min3)}");
Console.WriteLine();
Methods.DisplayMatrixLine(min4, "Метод наискорейшего спуска");
Console.WriteLine($"Минимальное значение: {MyFunction(min4)}");
Console.WriteLine();
Methods.DisplayMatrixLine(min5, "Метод сопряженных градиентов");
Console.WriteLine($"Минимальное значение: {MyFunction(min5)}");

// Console.ReadKey();
