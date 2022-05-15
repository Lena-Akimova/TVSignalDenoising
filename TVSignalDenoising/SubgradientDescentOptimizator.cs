using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

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
        /// <param name="_k">кол-во итераций</param>
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

        public (double[][], double[], int) Minimize()
        {
            Console.WriteLine($"\nСубградиентный спуск");
            var Xk = new double[k+1][];
            var Fk = new double[k+1];
            var Sk = new double[k+1][];
            var i = 0;
            Xk[i] = x0;
            Fk[i] = example.GetValueAt(Xk[i]).Value;
            var delta = eps + 1;
            try
            {
                while (delta > eps && i < k)
                {
                    Sk[i] = example.GetSubGradAt(Xk[i]);
                    //var SkNorm = Math.Sqrt(Sk[i].Sum(s => s * s));
                    var step = h / Math.Sqrt(i + 1);
                    //Console.WriteLine($"Шаг: {step}");
                    double[] toProj = Xk[i].Zip(Sk[i]).Select(xs => xs.First - step * xs.Second).ToArray();
                    Xk[i + 1] = Projection(toProj, example.BoxUp, example.BoxLow);
                    //Console.WriteLine($"X[{i + 1}] : {alglib.ap.format(Xk[i + 1], 3)}");
                    Fk[i + 1] = example.GetValueAt(Xk[i + 1]).Value;
                    //Console.WriteLine($"F[{i + 1}] : {Fk[i + 1]}\n ");
                    delta = Fk[i] - Fk[i + 1];
                    i++;
                }
            }
            catch (Exception ex)
            {
            }
            return (Xk, Fk, i);
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
