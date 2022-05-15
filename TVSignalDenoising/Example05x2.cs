using System;
using System.Collections.Generic;
using System.Text;

namespace TVSignalDenoising
{
    /// <summary>
    ///0.5x^2 + 5 -> min, 
    /// </summary>
    class Example05x2 : Example
    {
        /// <summary>
        /// <inheritdoc/>
        /// </summary>
        public override double[] GetSubGradAt(double[] x)
        {
            var s = new double[x.Length];
            for (int i = 0; i < s.Length; i++)
                s[i] = x[i];
            return s;
        }

        /// <summary>
        /// <inheritdoc/>
        /// </summary>
        public override double? GetValueAt(double[] x)
        {
            var fromDomain = x[0] <= BoxUp[0] && x[0] >= BoxLow[0];
            return fromDomain
                ? 0.5 * x[0] * x[0] + 5 
                : new double?();
        }

        /// <summary>
        /// <inheritdoc/>,
        ///  -||s||->min, x=[box],
        ///  -(Sum (si^2))^(1/2)->min, x=[box],
        ///  -(Sum (si^2))->min, x=[box]
        /// </summary>
        /// <returns></returns>
        public override double FindL()
        {
            double[,] A = new double[N, N]; //гессиан
            double[] B = new double[N]; // вектор коэфф-в b
            double[] scale = new double[N]; //масштаб

            for (int i = 0; i < N; i++)
            {
                A[i, i] = -2.0;
                scale[i] = 1; //один масштаб у всех иксов
            }

            double[] x;
            alglib.minqpstate state;
            alglib.minqpreport rep;
            alglib.minqpcreate(N, out state);

            alglib.minqpsetquadraticterm(state, A);
            alglib.minqpsetlinearterm(state, B);
            alglib.minqpsetscale(state, scale);


            alglib.minqpsetbc(state, BoxLow, BoxUp);
            alglib.minqpsetalgoquickqp(state, 0.0, 0.0, 0.0, 0, true);
            alglib.minqpoptimize(state);
            alglib.minqpresults(state, out x, out rep);

            var s = GetSubGradAt(x);
            var sum = 0.0;
            for (int i = 0; i < N; i++)
                sum += s[i] * s[i];
            sum = Math.Sqrt(sum);

            Console.WriteLine($" L = {sum} ");
            return sum;
        }

    }
}
