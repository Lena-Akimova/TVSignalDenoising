using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TVSignalDenoising
{
    /// <summary>
    /// 0.5*Sum(si-xi)^2 + lambda*Sum(|xi_1 - xi|) 
    /// </summary>
    class ExampleTV : Example
    {
        public double[] NoisedSignal { get; set; }
        public double[] DenoisedSignal { get; set; }
        public double[] Signal { get; set; }
        public double Lambda { get; set; }

        public ExampleTV()
        {
            Signal = FileIO.ReadSignal("", 50);
            NoisedSignal = NoiseSignal(Signal, 0.5, 0.05);
            N = Signal.Length;
        }

        /// <summary>
        /// Значение функции 0.5*Sum(si-xi)^2 + lambda*Sum(|xi_1 - xi|) в точке 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double? GetValueAt(double[] x)
        {
            var sumKvadr = NoisedSignal.Zip(x).Sum(sx => Math.Pow(sx.First - sx.Second, 2))*0.5;
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
                sub[i] = x[i] - NoisedSignal[i];
                if (i == 0)
                    sub[i] +=  - Lambda * Math.Sign(x[i + 1] - x[i]);
                else if (i == N - 1)
                    sub[i] += Lambda * Math.Sign(x[i] - x[i - 1]);
                else
                    sub[i] += - Lambda * Math.Sign(x[i + 1] - x[i]) + Lambda * Math.Sign(x[i] - x[i - 1]);
            }
            return sub;
        }


        /// <summary>
        /// <inheritdoc/>
        ///  (Sum (si^2))^(1/2)->max, x=[box],
        ///  (Sum (si^2))->max, x=[box]
        /// </summary>
        /// <returns></returns>
        public override double FindL()
        {
            ExampleSubgr subgr = new ExampleSubgr()
            {
                Lambda = this.Lambda,
                N = this.N,
                NoisedSignal = this.NoisedSignal,
                Scale = this.Scale,
                BoxLow = this.BoxLow,
                BoxUp=this.BoxUp
            };

            SubgradientDescentOptimizator so = new SubgradientDescentOptimizator(subgr, 100, NoisedSignal, 0.01, 0.0001);

            var res = so.Minimize();
            var MaxS = res.Item1[res.Item3];

            var sum = 0.0;
            for (int i = 0; i < N; i++)
                sum += MaxS[i] * MaxS[i];
            sum = Math.Sqrt(sum);

            Console.WriteLine($" L = {sum} ");
            return sum;
        }

        /// <summary>
        /// Среднеквадратическая ошибка
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
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

        public static double TV(double[] signal)
        {
            var tv = 0.0;
            for(int i=0; i<signal.Length-1; i++)
            {
                tv+=Math.Abs(signal[i + 1] - signal[i]);
            }
            return tv;
        }


    }
}
