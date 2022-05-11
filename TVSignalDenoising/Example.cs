using System;
using System.Collections.Generic;
using System.Text;

namespace TVSignalDenoising
{
    public abstract class Example
    {
        public double[] BoxUp { get; set; }
        public double[] BoxLow { get; set; }

        public int N { get; set; }

        public double Scale { get; set; }

        /// <summary>
        /// Субградиент в точке x
        public abstract double[] GetSubGradAt(double[] x);

        /// <summary>
        /// Значение функции в точке х
        /// </summary>
        public abstract double? GetValueAt(double[] x);

        /// <summary>
        /// ||s||-> max,
        /// x=[box]
        /// </summary>
        /// <returns></returns>
        public abstract double FindL();
    }
}
