using System;
using System.Collections.Generic;
using System.Text;

namespace TVSignalDenoising
{
    /// <summary>
    ///0.5(3*x1^2 + 10*x2^2) -> min, 
    /// </summary>
    class Example05x1x2 : Example
    {
        /// <summary>
        /// <inheritdoc/>
        /// </summary>
        public override double[] GetSubGradAt(double[] x)
        {
            var s = new double[x.Length];
            s[0] = 3 * x[0];
            s[1] = 10 * x[1];
            return s;
        }

        /// <summary>
        /// <inheritdoc/>
        /// </summary>
        public override double? GetValueAt(double[] x)
        {
            var fromDomain = true;
            return fromDomain
                ? 0.5*(3 * x[0] * x[0] + 10 * x[1] * x[1])
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
