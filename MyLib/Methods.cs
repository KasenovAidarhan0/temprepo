namespace MyLib;

public static class Methods
{
    private const int MaxIterations = 10000;
    
    public static Matrix ConjugateDescent(
        Func<Matrix, Matrix> fGrad,
        Matrix mQ,
        Matrix start,
        double eps
    )
    {
        Matrix x = start;
        Matrix h = (-1) * fGrad(x);
        for (int i = 0; i < MaxIterations; i++)
        {
            double beta = (fGrad(x) * mQ * h.Trsp()) /
                          (h * mQ * h.Trsp());
            
            double alpha = (-1) * (fGrad(x) * h.Trsp()) /
                           (h * mQ * h.Trsp());
            
            x += alpha * h;
            
            h = (-1) * fGrad(x) + beta * h;
            
            if (fGrad(x)[0, 0] < eps && fGrad(x)[0, 1] < eps)
            {
                break;
            }
        }
        return x;
    }
    
    public static double BisectionMethod(
        Func<double, double> f, 
        double a, 
        double b, 
        double eps)
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
    
    public static Matrix GradientDescent(
        Func<Matrix, Matrix> fGrad,
        Matrix mQ,
        Matrix start,
        double eps
    )
    {
        Matrix x = start;
        for (int i = 0; i < MaxIterations; i++)
        {
            double alpha = (fGrad(x) * fGrad(x).Trsp()) /
                           (fGrad(x) * mQ * fGrad(x).Trsp());
            x = x - alpha * fGrad(x);
            if (fGrad(x)[0, 0] < eps && fGrad(x)[0, 1] < eps)
            {
                break;
            }
        }
        return x;
    }
    
    public static Matrix GradientDescentDivide(
        Func<Matrix, double> f,
        Func<Matrix, Matrix> fGrad,
        Matrix start,
        double eps
    )
    {
        double alpha = 1;
        Matrix x = start;
        Matrix xPast = x;
        for (int i = 0; i < MaxIterations; i++)
        {

            for (int j = 0; j < MaxIterations; j++)
            {
                x = x - alpha * fGrad(x);
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
    
    public static Matrix GaussSeidel(
        Matrix mQ,
        Matrix mR,
        Matrix start,
        double eps
    )
    {
        double[,] e1 = { { 1, 0 } };
        double[,] e2 = { { 0, 1 } };
        double[][,] arr = [e1, e2];
        Matrix x = start;
        
        for (int i = 0; i < MaxIterations; i++)
        {
            var xPast = x;
            foreach (double[,] e in arr)
            {
                double alpha =  (-1) * (x * mQ * new Matrix(e).Trsp() + mR * new Matrix(e).Trsp()) /
                                (new Matrix(e) * mQ * new Matrix(e).Trsp());
                x = x + alpha * new Matrix(e);
            }

            if (EuclidNorm(x - xPast) < eps)
            {
                break;
            }
        }

        return x;
    }
    
    private static double EuclidNorm(Matrix a)
    {
        double res = 0;
        for (int i = 0; i < a.Rows; i++)
        {
            for (int j = 0; j < a.Columns; j++)
            {
                res += a[i, j] * a[i, j];
            }
        }

        return Math.Sqrt(res);
    }
    
    public static Matrix CoordinateDescent(
        Func<Matrix, double> f,
        Matrix start,
        double epsilon,
        double alpha
    )
    {
        double[,] e1 = { { 1, 0 } };
        double[,] e2 = { { 0, 1 } };
        Matrix x = start;
        for (int i = 0; i < MaxIterations; i++)
        {
            int conditionsFlag = 0;
            double[][,] eArr = [e1, e2];
            foreach (double[,] e in eArr)
            {
                Matrix stepXMinus = x - alpha * new Matrix(e);
                if (f(stepXMinus) < f(x))
                {
                    x = stepXMinus;
                }
                else
                {
                    conditionsFlag++;
                }

                Matrix stepXPlus = x + alpha * new Matrix(e);
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
            if (alpha < epsilon)
            {
                break;
            }
        }

        return x;
    }
}