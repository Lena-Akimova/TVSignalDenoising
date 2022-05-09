using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TVSignalDenoising
{
    class ExampleTV : Example
    {
        public override double[] BoxUp
        {
            get
            {
                var d = new double[N];
                for (int i = 0; i < N; i++)
                    d[i] = double.PositiveInfinity;
                return d;
            }
            set => BoxUp = value;
        }
        public override double[] BoxLow
        {
            get
            {
                var d = new double[N];
                for (int i = 0; i < N; i++)
                    d[i] = -double.PositiveInfinity;
                return d;
            }
            set => BoxLow = value;
        }
        public double[] NoisedSignal { get; set; }
        public override int N { get; set; }
        public double Lambda { get; set; }


        /// <summary>
        /// Значение функции Sum(si-xi)^2 + lambda*Sum(|xi_1 - xi|) в точке 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double? GetValueAt(double[] x)
        {
            var sumKvadr = NoisedSignal.Zip(x).Sum(sx => Math.Pow(sx.First - sx.Second, 2));
            var sumMod = 0.0;
            for (int i = 0; i < x.Length - 1; i++)
                sumMod += Math.Abs(x[i + 1] - x[i]);
            var inDomain = !x.Zip(BoxUp).Any(i => i.First > i.Second)&& !x.Zip(BoxLow).Any(i => i.First < i.Second);

            return inDomain
                ? sumKvadr + Lambda * sumMod
                : new double?();
        }

        /// <summary>
        /// Субградиент в точке
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double[] GetSubGradAt(double[] x)
        {
            var sub = new double[N];
            for (int i = 0; i < N; i++)
            {
                if (i == 0)
                {
                    var modPart = x[i + 1] == x[i] ? 0 : (x[i + 1] - x[i]) / Math.Abs(x[i + 1] - x[i]);
                    sub[i] = 2 * x[i] - 2 * NoisedSignal[i] - Lambda * modPart;
                }

                else if (i == N - 1)
                {
                    var modPart = x[i] == x[i - 1] ? 0 : (x[i] - x[i - 1]) / Math.Abs(x[i] - x[i - 1]);
                    sub[i] = 2 * x[i] - 2 * NoisedSignal[i] + Lambda * modPart;
                }
                else
                {
                    var modPart1 = x[i + 1] == x[i] ? 0 : (x[i + 1] - x[i]) / Math.Abs(x[i + 1] - x[i]);
                    var modPart2 = x[i] == x[i - 1] ? 0 : (x[i] - x[i - 1]) / Math.Abs(x[i] - x[i - 1]);
                    sub[i] = 2 * x[i] - 2 * NoisedSignal[i] + Lambda * modPart1 - Lambda * modPart2;
                }
            }

            return sub;
        }


        /// <summary>
        /// <inheritdoc/>
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
