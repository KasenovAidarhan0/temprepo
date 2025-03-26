namespace MyLib;

public struct Matrix
{
    private double[,] _data;
    public int Rows { get;  }
    public int Columns { get; }

    public Matrix(int rows, int columns)
    {
        Rows = rows;
        Columns = columns;
        _data = new double[rows, columns];
    }

    public Matrix(double[,] a)
    {
        Rows = a.GetUpperBound(0) + 1;
        Columns = a.GetUpperBound(1) + 1;
        _data = a;
    }

    public double[,] Data
    {
        get => _data;
        set
        {
            if (value.GetUpperBound(0) + 1 != Rows &&
                value.GetUpperBound(1) + 1 != Columns)
            {
                throw new Exception("Неподходящие размеры присваемового массива");
            }

            _data = value;
        }
    }

    public double this[int row, int col]
    {
        get => _data[row, col];
        set => _data[row, col] = value;
    }

    public void Print()
    {
        for (int i = 0; i < Rows; i++)
        {
            for (int j = 0; j < Columns; j++)
            {
                Console.Write("{0:f6}\t", _data[i, j]);
            }
            Console.WriteLine();
        }
    }

    public Matrix Trsp()
    {
        Matrix b = new Matrix(Columns, Rows);
        for (int i = 0; i < Rows; i++)
        {
            for (int j = 0; j < Columns; j++)
            {
                b[j, i] = _data[i, j];
            }
        }

        return b;
    }

    public int Count()
    {
        return Rows * Columns;
    }

    public static Matrix operator +(Matrix a, Matrix b)
    {
        if (a.Rows != b.Rows || a.Columns != b.Columns)
        {
            throw new Exception("Матрицы не одинакового размера.");
        }

        Matrix c = new Matrix(a.Rows, a.Columns);
        for (int i = 0; i < a.Rows; i++)
        {
            for (int j = 0; j < a.Columns; j++)
            {
                c[i, j] = a[i, j] + b[i, j];
            }
        }

        return c;
    }

    public static Matrix operator *(double a, Matrix b)
    {
        Matrix c = new Matrix(b.Rows, b.Columns);
        for (int i = 0; i < b.Rows; i++)
        {
            for (int j = 0; j < b.Columns; j++)
            {
                c[i, j] = a * b[i, j];
            }
        }

        return c;
    }

    public static Matrix operator *(Matrix a, Matrix b)
    {
        if (a.Columns != b.Rows)
        {
            throw new Exception("Матрицы не умножимы.");
        }

        Matrix c = new Matrix(a.Rows, b.Columns);
        for (int i = 0; i < a.Rows; i++)
        {
            for (int j = 0; j < b.Columns; j++)
            {
                c[i, j] = 0;
                for (int k = 0; k < a.Columns; k++)
                {
                    c[i, j] += a[i, k] * b[k, j];
                }
            }
        }

        return c;
    }

    public static Matrix operator -(Matrix a, Matrix b)
    {
        return a + (-1) * b;
    }

    public static double operator /(Matrix a, Matrix b)
    {
        if (a.Count() > 1 && a.Count() > 1)
        {
            throw new Exception("Матрицы должна иметь один элемент.");
        }
        return a[0, 0] / b[0, 0];
    }
}

