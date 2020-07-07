using System;
using System.Collections.Generic;
using System.Numerics;
using System.Linq;

// translation of ->  https://github.com/tyrex-team/benchmarks-attitude-smartphones/blob/master/src/Filters/Implementations/QMichelObsExtmagWtRep.m
namespace Logic.Ahrs.Algorithms.Tyrex.Implementations
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

        protected override Quaternion Update(Vector3 gyr, Vector3 acc, Vector3 mag, double dT)
        {
            var magUpdate = Math.Abs(mag.Length() - MagRefNorm) < MagNormThreshold;

            // % Do not consider a magUpdate for next [timeMagNopert] seconds if a ~magUpdate is detected
            LastMagPerturbation = MagUpdate ? LastMagPerturbation + dT : LastMagPerturbation;
            magUpdate = LastMagPerturbation < TimeMagNopert ? false : magUpdate;

            // % If we detect a perturbation, we replay all previous values from current time - [timeToSavePreviousData]
            // % to current time without magnetic field updates

            // obj.oldValues(end + 1,:) = [0 dT gyr acc mag obj.quaternion];
            // 1 + 1 + 3 + 3 + 3 + 4 -> 15 rows in oldValues
            if (MagUpdate && !magUpdate && OldValues.Count > 0)
            {
                var tmpBeta = Beta;
                Beta = 2;
                Quaternion = OldValues[0].Quaternion;

                OldValues.ForEach((value) =>
                {
                    var dTOld = value.DT;
                    var gyrOld = value.Gyr;
                    var accOld = value.Acc;
                    var magOld = value.Mag;

                    updateInterval(gyrOld, accOld, magOld, dTOld, magUpdate);
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
        void updateInterval(Vector3 gyr, Vector3 acc, Vector3 mag, double dT, bool magUpdate)
        {
            var q = Quaternion;

            acc /= acc.Length();
            mag /= mag.Length();

            var estimate_A = Quatrotate(q, AccRefNormalized);
            var estimate_M = Quatrotate(q, MagRefNormalized);
            
            var messure = CreateMatrix(acc, mag);
            var estimate = CreateMatrix(estimate_A, estimate_M);

            // delta = 2 * [obj.Ka * skew(estimate_A) ; obj.Km * magUpdate * skew(estimate_M)]';
            var a = Skew(estimate_A).Multiply(Ka);
            Matrix4x4? m = null;
            if(magUpdate)
            {
                m = Skew(estimate_M).Multiply(Km);
            }

            Matrix4x4 am = a; // []
            if (m.HasValue)
            {
                am = Matrix4x4.Add(a, m.Value);
            }
            
            var delta = am.Transpose() * 2;
            //

            var eye3 = Matrix4x4.Identity; // is diagonalidentity

            var dq = (messure - estimate) * ((delta * delta.Transpose() * eye3.Add(0.00001)).Inverse() * delta).Transpose();

            var gyrQ = new Quaternion(gyr.X, gyr.Y, gyr.Z, 0);
            var dqAsQ = System.Numerics.Quaternion.CreateFromRotationMatrix(dq); // matlab code seems to assume as if dq is a vector

            // does order of (Quaternin)muliplication matter ..?
            //var qDot = 0.5 * q * gyroQ + Beta + q * dqAsQ;
            var qDot = (q * gyrQ).Multiply(0.5f) + (q * dqAsQ).Multiply(Beta);
            q += qDot.Multiply(dT);
            q = q.Divide(q.Length());

            Quaternion = q;
            Logic.Utils.Log.Message("QMichelObsExtmagWtRep quaternion became: " + Quaternion);
        }

    }
    public class Value
    {
        public double LastMagPerturbation { get; set; }
        public double DT { get; set; }
        public Vector3 Gyr { get; set; }
        public Vector3 Acc { get; set; }
        public Vector3 Mag { get; set; }
        public Quaternion Quaternion { get; set; }

        public Value(double lastMagPerturbation, double dT, Vector3 gyr, Vector3 acc, Vector3 mag, Quaternion quaternion)
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
