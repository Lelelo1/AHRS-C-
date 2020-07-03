using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Spatial.Euclidean;
// https://github.com/tyrex-team/benchmarks-attitude-smartphones/blob/master/src/Filters/Implementations/QMichelObsExtmagWtRep.m
namespace Logic.Ahrs.Magnetometer.Tyrex.Implementations
{
    public class QMichelObsExtmagWtRep : AttitudeFilter
    {

        public 

        double MagNormThreshold { get; set; } = 15;
        double Beta { get; set; } = 0.3;
        double Ka { get; set; } = 2;
        double Km { get; set; } = 1;

        // oldValues = [];
        double LastMagPerturbation { get; set; }
        // magUpdate;

        double TimeMagNopert { get; set; } = 2;
        double TimeToSavePreviousData { get; set; } = 3;

        public override Quaternion Update(Vector<double> gyr, Vector<double> acc, Vector<double> mag, double dT)
        {
            
        }
    }
}
