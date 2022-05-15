using System;
using System.Collections.Generic;
using System.Text;

namespace TVSignalDenoising
{
    interface IOptimizator
    {
        /// <summary>
        /// Минимизация
        /// </summary>
        /// <returns>Кортеж из массива точек и значений функции, итерации с ответом</returns>
        (double[][], double[], int) Minimize();
    }
}
