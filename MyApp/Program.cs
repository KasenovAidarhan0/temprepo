using MyLib;

// Функция
double MyFunction(Matrix a)
{
    double x1 = a[0, 0];
    double x2 = a[0, 1];
    return 4 * x1 * x1 + x2 * x2 - x1 * x2 + 7 * x1 + 3 * x2;
}

// Градиент
Matrix MyGrad(Matrix a)
{
    double x1 = a[0, 0];
    double x2 = a[0, 1];
    Matrix res = new Matrix(1, 2)
    {
        [0, 0] = 8 * x1 - x2 + 7,
        [0, 1] = 2 * x2 - x1 + 3
    };
    return res;
}

// Q
var mQ = new Matrix(2, 2)
{
    Data = new double[,]
    {
        { 8, -1 },
        { -1, 2 }
    }
};

// R
var mR = new Matrix(1, 2)
{
    Data = new double[,]
    {
        { 7, 3 }
    }
};

// Точность
double tolerance = 0.001;

// Нaчальная точка
Matrix start = new Matrix(1, 2) { Data = new double[,] { { 0, 0 } } };

// Метод координатного спуска 
var min1 = Methods.CoordinateDescent(MyFunction, start, tolerance, 1);
// Метод Гаусса-Зейделя
var min2 = Methods.GaussSeidel(mQ, mR, start, tolerance);
// Метод градиентного спуска
var min3 = Methods.GradientDescentDivide(MyFunction, MyGrad, start, tolerance);
// Метод наискорейшего спуска
var min4 = Methods.GradientDescent(MyGrad, mQ, start, tolerance);
// Метод сопряженных градиентов
var min5 = Methods.ConjugateDescent(MyGrad, mQ, start, tolerance);

min1.Print();
min2.Print();
min3.Print();
min4.Print();
min5.Print();

// Console.ReadKey();
