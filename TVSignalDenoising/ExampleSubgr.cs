using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TVSignalDenoising
{
    /// <summary>
    /// - Скалярный квадрат субградиента ф-ии с l1-регуляризацией
    /// </summary>
    class ExampleSubgr : ExampleTV
    {
        public override double FindL()
        {
            throw new NotImplementedException();
        }

        public override double[] GetSubGradAt(double[] x)
        {
            var sub = new double[N];
            for (int i = 0; i < N; i++)
            {
                if (i == 0)
                {
                    sub[i] = -(2 * (x[i] - NoisedSignal[i] - Lambda * Math.Sign(x[i + 1] - x[i]))
                        * (1 + Lambda * (1 / Math.Abs(x[i + 1] - x[i]) - Math.Sign(x[i + 1] - x[i]) / (x[i + 1] - x[i])))
                        -2 * Lambda * (x[i + 1] - NoisedSignal[i + 1] - Lambda * Math.Sign(x[i + 2] - x[i + 1]) + Lambda * Math.Sign(x[i + 1] - x[i]))
                        * (1 / Math.Abs(x[i + 1] - x[i]) - Math.Sign(x[i + 1] - x[i]) / (x[i + 1] - x[i])));
                }
                else if (i == N - 1)
                {
                    sub[i] = -(2 * (x[i] - NoisedSignal[i] + Lambda * Math.Sign(x[i] - x[i - 1]))
                        * (1 + Lambda * (1 / Math.Abs(x[i] - x[i - 1]) - Math.Sign(x[i] - x[i - 1]) / (x[i] - x[i - 1])))
                        -2 * Lambda * (x[i - 1] - NoisedSignal[i - 1] - Lambda * Math.Sign(x[i] - x[i - 1]) + Lambda * Math.Sign(x[i - 1] - x[i - 2]))
                        * (1 / Math.Abs(x[i] - x[i - 1]) - Math.Sign(x[i] - x[i - 1]) / (x[i] - x[i - 1])));
                }
                else
                {
                    sub[i] = -(-2 * (x[i - 1] - NoisedSignal[i - 1] - Lambda * Math.Sign(x[i] - x[i - 1]))
                        * Lambda * (1 / Math.Abs(x[i] - x[i - 1]) - Math.Sign(x[i] - x[i - 1]) / (x[i] - x[i - 1]))
                        + 2 * (x[i] - NoisedSignal[i] - Lambda * Math.Sign(x[i + 1] - x[i]) + Lambda * Math.Sign(x[i] - x[i - 1]))

                        * (1 + Lambda * (1 / Math.Abs(x[i + 1] - x[i]) - Math.Sign(x[i + 1] - x[i]) / (x[i + 1] - x[i])))
                        + Lambda * (1 / Math.Abs(x[i] - x[i - 1]) - Math.Sign(x[i] - x[i - 1]) / (x[i] - x[i - 1]))
                        - 2 * (x[i + 1] - NoisedSignal[i + 1] + Lambda * Math.Sign(x[i + 1] - x[i]))
                        * Lambda * (1 / Math.Abs(x[i + 1] - x[i]) - Math.Sign(x[i + 1] - x[i]) / (x[i + 1] - x[i])));
                }
            }
            return sub;
        }


        /// <summary>
        /// - Скалярный квадрат субградиента
        /// -((x1-s1-lambda*sign(x2-x1))^2 +
        /// + (x2-s2-lambda*sign(x3-x2) + lambda*sign(x2-x1))^2 +...
        /// + (xi-si-lambda*sign(xi1-xi) + lambda*sign(xi-xi_1))^2 +...
        /// + (xi_n-si_n-lambda*sign(xi_n-xi_n_1))^2)
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double? GetValueAt(double[] x)
        {
            var sub = ((ExampleTV)this).GetSubGradAt(x);
            return sub.Sum(s => -s * s);
        }
    }
}
