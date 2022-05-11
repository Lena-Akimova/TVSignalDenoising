using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TVSignalDenoising
{
    /// <summary>
    /// Sum(si-xi)^2 + lambda*Sum(|xi_1 - xi|) 
    /// </summary>
    class ExampleTV : Example
    {
        public double[] NoisedSignal { get; set; }
        public double[] DenoisedSignal { get; set; }
        public double[] Signal { get; set; }
        public double Lambda { get; set; }

        public ExampleTV()
        {
            Signal = FileIO.ReadSignal("", 128);
            NoisedSignal = NoiseSignal(Signal, 0.5, 0.05);
            N = Signal.Length;
        }

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
                    var modPart = x[i + 1] == x[i] ? 0 : x[i + 1] > x[i] ? 1 : -1; //(x[i + 1] - x[i]) / Math.Abs(x[i + 1] - x[i]);
                    sub[i] = 2 * x[i] - 2 * NoisedSignal[i] - Lambda * modPart;
                }

                else if (i == N - 1)
                {
                    var modPart = x[i] == x[i - 1] ? 0 : x[i] > x[i - 1] ? 1 : -1;// (x[i] - x[i - 1]) / Math.Abs(x[i] - x[i - 1]);
                    sub[i] = 2 * x[i] - 2 * NoisedSignal[i] + Lambda * modPart;
                }
                else
                {
                    var modPart1 = x[i + 1] == x[i] ? 0 : x[i + 1] > x[i] ? 1 : -1;// (x[i + 1] - x[i]) / Math.Abs(x[i + 1] - x[i]);
                    var modPart2 = x[i] == x[i - 1] ? 0 : x[i] > x[i - 1] ? 1 : -1;// (x[i] - x[i - 1]) / Math.Abs(x[i] - x[i - 1]);
                    sub[i] = 2 * x[i] - 2 * NoisedSignal[i] - Lambda * modPart1 + Lambda * modPart2;
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

        public static double RMSE(double[] x, double[] y)
        {
            var s = x.Zip(y).Sum(xy => Math.Pow(xy.First - xy.Second, 2));
            return Math.Sqrt(s / x.Length);
        }

        public static double[] NoiseSignal(double[] signal, double k = 0.3, double prob = 0.03)
        {
            var noised = new double[signal.Length];
            Random r = new Random();
            for (var i = 0; i < signal.Length; i++)
            {
                var rand = (r.NextDouble() * 2 - 1) * k;
                if (r.NextDouble() < prob)
                    noised[i] = signal[i] + rand * 7;
                else
                    noised[i] = signal[i] + rand;
            }
            return noised;
        }


    }
}
