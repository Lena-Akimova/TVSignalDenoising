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
        /// <returns>Кортеж из массива точек и значений функции, итерации с ответом, кол-ва шагов</returns>
        (double[][], double[], int, int) Minimize();
    }
}
