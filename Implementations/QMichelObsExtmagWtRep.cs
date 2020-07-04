using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Spatial.Euclidean;


//using MathNet.Spatial.Euclidean;
// https://github.com/tyrex-team/benchmarks-attitude-smartphones/blob/master/src/Filters/Implementations/QMichelObsExtmagWtRep.m
namespace Logic.Ahrs.Magnetometer.Tyrex.Implementations
{
    public class QMichelObsExtmagWtRep : AttitudeFilter
    {

        double MagNormThreshold { get; set; } = 15;

        double Beta { get; set; } = 0.3;
        double Ka { get; set; } = 2;
        double Km { get; set; } = 1;



        List<Value> OldValues { get; set; } = new List<Value>();
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

            // obj.oldValues(end + 1,:) = [0 dT gyr acc mag obj.quaternion];
            // 1 + 1 + 3 + 3 + 3 + 4 -> 15 rows in oldValues
            if (MagUpdate && !magUpdate)
            {
                var tmpBeta = Beta;
                Beta = 2;
                Quaternion = OldValues[0].Quaternion;

                OldValues.ForEach((value) =>
                {
                    var dT = value.DT;
                    var gyrOld = value.Gyr;
                    var accOld = value.Acc;
                    var magOld = value.Mag;

                    updateInterval(gyrOld, accOld, magOld, dT, magUpdate);
                });
                Beta = tmpBeta;
            }

            // % Save values in case of reprocessing
            OldValues.Add(new Value(0, dT, gyr, acc, mag, Quaternion));
            OldValues.ForEach((value) =>
            {
                value.LastMagPerturbation += dT;

            });
            OldValues.RemoveAll((v) => v.LastMagPerturbation > TimeToSavePreviousData);

            // % Use the normal filter
            updateInterval(gyr, acc, mag, dT, magUpdate);
            MagUpdate = magUpdate;

            return Quaternion;
        }
        void updateInterval(Vector<double> gyr, Vector<double> acc, Vector<double> mag, double dT, bool magUpdate)
        {
            var q = Quaternion;

            acc /= acc.L2Norm();
            mag /= mag.L2Norm();

            var estimate_A = Quatrotate(q, AccRefNormalized);
            var estimate_M = Quatrotate(q, MagRefNormalized);

            var messure = Matrix<double>.Build.DenseOfRowVectors(acc, mag);
            var estimate = Matrix<double>.Build.DenseOfRowVectors(estimate_A, estimate_M);

            var a = Ka * Skew(estimate_A);
            var m = magUpdate ? Km * Skew(estimate_M) : null; ;

            Matrix<double> delta = a;
            if (m != null)
            {
                delta.Add(m);
            }
            delta *= 2;

            var eye3 = Matrix<double>.Build.DiagonalIdentity(3, 3);
            var dq = (messure - estimate) * ((delta * delta.Transpose() + 100000 * eye3).Power(-1) * delta).Transpose();

            var gyroQ = new Quaternion(0, gyr[0], gyr[1], gyr[2]);
            var dqAsQ = new Quaternion(0, dq.Row(0)[0], dq.Row(0)[1], dq.Row(0)[2]); // should be a be a vector;

            var qDot = 0.5 * q * gyroQ + Beta + q * dqAsQ;
            q += qDot * dT;
            q /= q.Norm;

            Quaternion = q;

        }
    }
    public class Value
    {
        public double LastMagPerturbation { get; set; }
        public double DT { get; set; }
        public Vector<double> Gyr { get; set; }
        public Vector<double> Acc { get; set; }
        public Vector<double> Mag { get; set; }
        public Quaternion Quaternion { get; set; }

        public Value(double lastMagPerturbation, double dT, Vector<double> gyr, Vector<double> acc, Vector<double> mag, Quaternion quaternion)
        {
            LastMagPerturbation = lastMagPerturbation;
            DT = dT;
            Gyr = gyr;
            Acc = acc;
            Mag = mag;
            Quaternion = quaternion;
        }
    }
}
