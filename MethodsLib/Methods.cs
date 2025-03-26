namespace MethodsLib;

public static class Methods
{
    public static double[,] ConjugateDescent(
        Func<double[,], double[,]> fGrad,
        double[,] mQ,
        double[,] start,
        double eps
    )
    {
        double[,] x = start;
        double[,] h = MatrixScalarMultiplication(fGrad(x), -1);
        int maxIterations = 10000;
        for (int i = 0; i < maxIterations; i++)
        {
            double beta = Multiply(Multiply(fGrad(x), mQ), Transpose(h))[0, 0] /
                          Multiply(Multiply(h, mQ), Transpose(h))[0, 0];
            h = MatrixAddition(MatrixScalarMultiplication(fGrad(x), -1), MatrixScalarMultiplication(h, beta));
            
            double[,] alphaUp = Multiply(fGrad(x), Transpose(h));
            double[,] alphaDown = Multiply(Multiply(h, mQ), Transpose(h));
            double alpha = (-1) * alphaUp[0, 0] / alphaDown[0, 0];
            
            x = MatrixAddition(x, MatrixScalarMultiplication(h, alpha));
            if (fGrad(x)[0, 0] < eps && fGrad(x)[0, 1] < eps)
            {
                break;
            }
        }
        return x;
    }
    
    public static double[,] GradientDescent(
        Func<double[,], double[,]> fGrad,
        double[,] mQ,
        double[,] start,
        double eps
    )
    {
        double[,] x = start;
        int maxIterations = 10000;
        for (int i = 0; i < maxIterations; i++)
        {
            double a = Multiply(fGrad(x), Transpose(fGrad(x)))[0, 0];
            double b = Multiply(Multiply(fGrad(x), mQ), Transpose(fGrad(x)))[0, 0];
            double alpha = a / b;
            x = MatrixSubtraction(x, MatrixScalarMultiplication(fGrad(x), alpha));
            if (fGrad(x)[0, 0] < eps && fGrad(x)[0, 1] < eps)
            {
                break;
            }
        }
        return x;
    }
    
    public static double[,] GradientDescentDivide(
        Func<double[,], double> f,
        Func<double[,], double[,]> fGrad,
        double[,] start,
        double eps
    )
    {
        int maxIterations = 10000;
        double alpha = 1;
        double[,] x = start;
        double[,] xPast = x;
        for (int i = 0; i < maxIterations; i++)
        {

            for (int j = 0; j < maxIterations; j++)
            {
                x = MatrixSubtraction(x, MatrixScalarMultiplication(fGrad(x), alpha));
                if (f(x) < f(xPast))
                {
                    break;
                }
                x = xPast;
                alpha /= 2;
            }
            xPast = x;
            if (fGrad(x)[0, 0] < eps && fGrad(x)[0, 1] < eps)
            {
                break;
            }
        }

        return x;
    }
    
    public static double[,] GaussSeidel(
        double[,] mQ,
        double[,] mR,
        double[,] start,
        double eps
    )
    {
        double[,] e1 = { { 1, 0 } };
        double[,] e2 = { { 0, 1 } };
        double[][,] arr = [e1, e2];
        double[,] x = start;

        int maxIterations = 10000;
        for (int i = 0; i < maxIterations; i++)
        {
            var xPast = x;
            foreach (double[,] e in arr)
            {
                double[,] eT = Transpose(e);
                double a = Multiply(Multiply(x, mQ), eT)[0, 0] + Multiply(mR, eT)[0, 0];
                double b = Multiply(Multiply(e, mQ), eT)[0, 0];
                double alpha = (-1) * a / b;
                x = MatrixAddition(x, MatrixScalarMultiplication(e, alpha));
            }

            if (EuclidNorm(MatrixSubtraction(x, xPast)) < eps)
            {
                break;
            }
        }

        return x;
    }

    private static double EuclidNorm(double[,] a)
    {
        int rows = a.GetUpperBound(0) + 1;
        int columns = a.GetUpperBound(1) + 1;
        double res = 0;
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                res += a[i, j] * a[i, j];
            }
        }

        return Math.Sqrt(res);
    }

    public static double BisectionMethod(Func<double, double> f, double a, double b, double eps)
    {
        while (Math.Abs(b - a) > eps)
        {
            double mid = (a + b) / 2;
            double delta = (b - a) / 4;
            if (f(a + delta) < f(mid))
            {
                b = mid;
            }
            else if (f(b - delta) < f(mid))
            {
                a = mid;
            }
            else
            {
                a += delta;
                b -= delta;
            }
        }

        return (a + b) / 2;
    }
    
    public static double[,] CoordinateDescent(
        Func<double[,], double> f, 
        double[,] start, 
        double eps,
        double alpha)
    {
        int maxIterations = 10000;
        double[,] e1 = { { 1, 0 } };
        double[,] e2 = { { 0, 1 } };
        double[][,] arr = [e1, e2];
        double[,] x = start;

        for (int i = 0; i < maxIterations; i++)
        {
            int conditionsFlag = 0;
            foreach (double[,] e in arr)
            {
                double[,] stepXMinus = MatrixSubtraction(x, MatrixScalarMultiplication(e, alpha));
                if (f(stepXMinus) < f(x))
                {
                    x = stepXMinus;
                }
                else
                {
                    conditionsFlag++;
                }

                double[,] stepXPlus = MatrixAddition(x, MatrixScalarMultiplication(e, alpha));
                if (f(stepXPlus) < f(x))
                {
                    x = stepXPlus;
                }
                else
                {
                    conditionsFlag++;
                }
            }

            if (conditionsFlag == 4)
            {
                alpha /= 2;
            }

            if (alpha < eps)
            {
                break;
            }
        }

        return x;
    }

    public static double[,] MatrixScalarMultiplication(double[,] matrix, double c)
    {
        int rows = matrix.GetUpperBound(0) + 1;
        int columns = matrix.GetUpperBound(1) + 1;
        double[,] result = new double[rows, columns];

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                result[i, j] = matrix[i, j] * c;
            }
        }

        return result;
    }

    public static double[,] MatrixAddition(double[,] matrix1, double[,] matrix2)
    {
        int rows1 = matrix1.GetUpperBound(0) + 1;
        int columns1 = matrix1.GetUpperBound(1) + 1;
        
        int rows2 = matrix2.GetUpperBound(0) + 1;
        int columns2 = matrix2.GetUpperBound(1) + 1;

        if (rows1 != rows2 || columns1 != columns2)
        {
            throw new Exception("Матрицы не имеют одинакового размера");
        }

        double[,] result = new double[rows1, columns1];

        for (int i = 0; i < rows1; i++)
        {
            for (int j = 0; j < columns1; j++)
            {
                result[i, j] = matrix1[i, j] + matrix2[i, j];
            }
        }

        return result;
    }

    public static double[,] MatrixSubtraction(double[,] matrix1, double[,] matrix2)
    {
        double[,] invMatrix = MatrixScalarMultiplication(matrix2, -1);
        return MatrixAddition(matrix1, invMatrix);
    }

    public static void DisplayMatrix(double[,] matrix, string name)
    {
        int rows = matrix.GetUpperBound(0) + 1;
        int columns = matrix.GetUpperBound(1) + 1;
        Console.WriteLine($"matrix: {name}");
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                Console.Write($"{matrix[i, j]}\t");
            }
            Console.WriteLine();
        }
        Console.WriteLine();
    }

    public static void DisplayMatrixLine(double[,] matrix, string name)
    {
        if (matrix.GetUpperBound(0) + 1 > 1)
        {
            throw new Exception("Матрица имеет больше одной строки.");
        }
        Console.Write($"{name}: ");
        for (int i = 0; i < matrix.GetUpperBound(1) + 1; i++)
        {
            Console.Write($"{matrix[0, i]}\t");
        }
        Console.WriteLine();
    }

    public static double[,] Transpose(double[,] matrix)
    {
        int rows = matrix.GetUpperBound(0) + 1;
        int columns = matrix.GetUpperBound(1) + 1;
        double[,] result = new double[columns, rows];
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                result[j, i] = matrix[i, j];
            }
        }

        return result;
    }

    public static double[,] Multiply(double[,] a, double[,] b)
    {
        int aRows = a.GetUpperBound(0) + 1;
        int aColumns = a.GetUpperBound(1) + 1;

        int bRows = b.GetUpperBound(0) + 1;
        int bColumns = b.GetUpperBound(1) + 1;

        if (aColumns != bRows)
        {
            throw new Exception("Матрицы не перемножаются");
        }

        double[,] result = new double[aRows, bColumns];

        for (int i = 0; i < aRows; i++)
        {
            for (int j = 0; j < bColumns; j++)
            {
                result[i, j] = 0;
                for (int k = 0; k < aColumns; k++)
                {
                    result[i, j] += a[i, k] * b[k, j];
                }
            }
        }

        return result;
    }
}