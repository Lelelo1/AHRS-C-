using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Spatial.Euclidean;
// https://github.com/tyrex-team/benchmarks-attitude-smartphones/blob/master/src/Filters/Implementations/QMichelObsExtmagWtRep.m
namespace Logic.Ahrs.Magnetometer.Tyrex.Implementations
{
    public class QMichelObsExtmagWtRep : AttitudeFilter
    {

        double MagNormThreshold { get; set; } = 15;
        double Beta { get; set; } = 0.3;
        double Ka { get; set; } = 2;
        double Km { get; set; } = 1;

        // oldValues = [];
        double LastMagPerturbation { get; set; }
        bool MagUpdate { get; set; }

        double TimeMagNopert { get; set; } = 2;
        double TimeToSavePreviousData { get; set; } = 3;

        public QMichelObsExtmagWtRep()
        {
            LastMagPerturbation = TimeMagNopert;
            MagUpdate = false;
        }

        public override Quaternion Update(Vector<double> gyr, Vector<double> acc, Vector<double> mag, double dT)
        {
            var magUpdate = Math.Abs(mag.L2Norm() - MagRefNorm) < MagNormThreshold;


            // % Do not consider a magUpdate for next [timeMagNopert] seconds if a ~magUpdate is detected
            LastMagPerturbation = MagUpdate ? LastMagPerturbation + dT : LastMagPerturbation;
            magUpdate = LastMagPerturbation < TimeMagNopert ? false : magUpdate;

            // % If we detect a perturbation, we replay all previous values from current time - [timeToSavePreviousData]
            // % to current time without magnetic field updates

            
            if (MagUpdate && !magUpdate)
            {
                var tmpBeta = Beta;
                Beta = 2;
                Quaternion = 
            }

        }
    }
}
