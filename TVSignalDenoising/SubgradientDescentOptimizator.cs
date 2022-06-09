using System;
using System.Linq;

namespace TVSignalDenoising
{
    /// <summary>
    /// Субградиентный метод
    /// </summary>
    class SubgradientDescentOptimizator : IOptimizator
    {
        private Example example;
        private double[] x0;
        private double h;
        private int k;
        private double eps;

        /// <param name="ex">функция</param>
        /// <param name="_k">количество итераций</param>
        /// <param name="_x0">начальная точка</param>
        /// <param name="_h">начальный шаг</param>
        public SubgradientDescentOptimizator(Example ex, int _k, double[] _x0, double _h, double _eps)
        {
            example = ex;
            x0 = _x0;
            h = _h;
            eps = _eps;
            k = _k;
        }

        public (double[][], double[], int, int) Minimize()
        {
            var Xk = new double[k+1][];
            var Fk = new double[k+1];
            var Sk = new double[k+1][];
            var i = 0;
            Xk[i] = x0;
            Fk[i] = example.GetValueAt(Xk[i]).Value;
            var delta = eps + 1;
            while (true)
            {
                try
                {
                    while (Math.Abs(delta )> eps && i < k)
                    {
                        Sk[i] = example.GetSubGradAt(Xk[i]);
                        var step = h / Math.Sqrt(i + 1);
                        double[] toProj = Xk[i].Zip(Sk[i]).Select(xs => xs.First - step * xs.Second).ToArray();
                        Xk[i + 1] = Projection(toProj, example.BoxUp, example.BoxLow);
                        Fk[i + 1] = example.GetValueAt(Xk[i + 1]).Value;
                        delta = Fk[i] - Fk[i + 1];
                        i++;
                    }
                    if (delta > eps)
                    {
                        var Xk1 = new double[2 * k + 2][];
                        var Fk1 = new double[2 * k + 2];
                        var Sk1 = new double[2 * k + 2][];
                        Array.Copy(Xk, Xk1, Xk.Length);
                        Array.Copy(Fk, Fk1, Fk.Length);
                        Array.Copy(Sk, Sk1, Sk.Length);
                        Xk = Xk1;
                        Fk = Fk1;
                        Sk = Sk1;
                        k=2*k;
                    }
                    else break;
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                    break;
                }
                
            }

            return (Xk, Fk, i,i+1);
        }


        /// <summary>
        /// Проекция точки на множество
        /// Задача условной квадратичной минимизации
        /// ||x-x0||^2 -> min = Sum((xi-x0i)^2)-> min
        /// xi>=0
        /// </summary>
        /// <param name="x0">проектируемая точка</param>
        /// <returns></returns>
        public double[] Projection(double[] x0, double[] boxUp, double[] boxLow)
        {
            int n = example.N;

            double[,] A = new double[n, n];
            double[] B = new double[n];
            double[] scale = new double[n];

            for (int i = 0; i < n; i++)
            {
                A[i, i] = 2.0; 
                B[i] = -2 * x0[i];
                scale[i] = example.Scale; 
            }

            double[] x;
            alglib.minqpstate state;
            alglib.minqpreport rep;

            alglib.minqpcreate(n, out state);
            alglib.minqpsetquadraticterm(state, A);
            alglib.minqpsetlinearterm(state, B);
            alglib.minqpsetbc(state, boxLow, boxUp);
            alglib.minqpsetscale(state, scale);
            alglib.minqpsetalgoquickqp(state, 0.0, 0.0, 0.0, 0, true);
            alglib.minqpoptimize(state);
            alglib.minqpresults(state, out x, out rep);
            return x;
        }

    }
}
